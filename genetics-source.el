;;; genetics-source.el --- Swappable kit sources: native parser or genome-cli -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Package-Requires: ((emacs "29.1"))
;; Keywords: data, tools
;; URL: https://gitlab.com/davidawad/genetics-el

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; A kit can come from two interchangeable sources, chosen by
;; `genetics-source-function':
;;
;; - `genetics-source-native': the Emacs Lisp parser in genetics-parse.el
;;   (array files and VCFs).
;; - `genetics-source-genome-cli': the local `genome' executable
;;   (genome-cli), which parses and normalizes every format and answers
;;   with genome/v1 JSON envelopes.  Its kits, summaries and genotype
;;   records are mapped onto `genetics-kit' and `genetics-snp', so
;;   browse, lookup, report, compare and export work on either source.
;;
;; The default, `genetics-source-auto', uses genome-cli when
;; `genetics-genome-executable' is found and the native parser otherwise.
;;
;; genome-cli commands used (all with --format json, stdout only):
;;
;;   genome import FILE                       kind "kits"
;;   genome summary ID                        kind "summary"
;;   genome lookup ID --rsid RSID             kind "genotypes"
;;   genome export ID --region C[:S-E]
;;                --limit N --offset M        kind "genotypes"
;;   genome compare ID-A ID-B                 kind "compare"
;;
;; Every command has a pure twin that returns the argv without running
;; it: `genetics-genome-argv', `genetics-source-genome-cli-explain'.
;; genome-cli runs locally; this file opens no network connection.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)
(require 'genetics-core)
(require 'genetics-parse)
(require 'genetics-stats)
(require 'genetics-annotate)

;;;; Customization

(defcustom genetics-genome-executable "genome"
  "Name or absolute path of the genome-cli executable."
  :type 'string
  :group 'genetics)

(defcustom genetics-source-function #'genetics-source-auto
  "Function that turns a file into a `genetics-kit'.
It is called with the file name and the keyword arguments of
`genetics-parse-file'.  `genetics-source-auto' uses genome-cli when its
executable is available and the native parser otherwise."
  :type '(choice (const :tag "genome-cli when installed, else native"
                        genetics-source-auto)
                 (const :tag "genome-cli" genetics-source-genome-cli)
                 (const :tag "Native Emacs Lisp parser" genetics-source-native)
                 function)
  :group 'genetics)

(defcustom genetics-genome-page-size 5000
  "Records requested per `genome query' call when walking a kit."
  :type 'integer
  :group 'genetics)

(defconst genetics-genome-schema "genome/v1"
  "The genome-cli JSON contract this file understands.")

;;;; Choosing a source

(defun genetics-genome-available-p ()
  "Return the genome-cli executable path, or nil if it is not installed."
  (executable-find genetics-genome-executable))

(defun genetics-source-native (file &rest args)
  "Parse FILE in Emacs Lisp; ARGS as for `genetics-parse-file'."
  (apply #'genetics-parse-file file args))

(defun genetics-source-auto (file &rest args)
  "Open FILE with genome-cli if it is installed, else natively.
ARGS are passed on unchanged."
  (apply (if (genetics-genome-available-p)
             #'genetics-source-genome-cli
           #'genetics-source-native)
         file args))

(defun genetics-source-open (file &rest args)
  "Return a kit for FILE from `genetics-source-function'.
ARGS are the keyword arguments of `genetics-parse-file'."
  (apply genetics-source-function file args))

;;;; argv (pure)

(defun genetics-genome-argv (command &rest args)
  "Return the full genome-cli argv for COMMAND with ARGS.
COMMAND is a symbol: `import' (FILE), `summary' (ID), `lookup' (ID
RSID), `query' (ID CHROM OFFSET &optional START END) or `compare' (ID-A
ID-B).  The first element is `genetics-genome-executable'.  Pure: nothing
is run."
  (cons genetics-genome-executable
        (append
         (pcase command
           ('import (append (list "import" (expand-file-name (car args)))
                            (when (cadr args) (list "--replace"))))
           ('kits (list "kits"))
           ('summary (list "summary" (car args)))
           ('lookup (list "lookup" (car args) "--rsid" (cadr args)))
           ('query
            ;; genome-cli pages a region through `export' (kind "genotypes").
            (pcase-let ((`(,id ,chrom ,offset ,start ,end) args))
              (list "export" id
                    "--region" (if start (format "%s:%d-%d" chrom start end) chrom)
                    "--limit" (number-to-string genetics-genome-page-size)
                    "--offset" (number-to-string offset))))
           ('compare (list "compare" (car args) (cadr args)))
           (_ (genetics--error 'genetics-genome-error
                               "Unknown genome command %S" command)))
         (list "--format" "json"))))

(defun genetics-genome-argv-string (argv)
  "Return ARGV as one shell-quoted command line."
  (mapconcat #'shell-quote-argument argv " "))

(defun genetics-source-genome-cli-explain (file)
  "Return the exact command `genetics-source-genome-cli' would run for FILE.
Pure twin of `genetics-source-genome-cli': nothing is executed.  The
import is followed by `genome summary ID --format json' for the kit
id it returns.  Interactively, show the command in the echo area."
  (interactive (list (read-file-name "Explain genome-cli import of: "
                                     genetics-data-directory nil t)))
  (let ((cmd (genetics-genome-argv-string (genetics-genome-argv 'import file))))
    (when (called-interactively-p 'interactive)
      (message "%s" cmd))
    cmd))

;;;; Running genome-cli

(defun genetics--genome-parse-envelope (text argv)
  "Parse genome/v1 envelope TEXT printed by ARGV; return the alist.
Signal `genetics-genome-error' when TEXT is not a genome/v1 envelope."
  (let ((env (condition-case nil
                 (json-parse-string (string-trim text)
                                    :object-type 'alist :array-type 'list
                                    :null-object nil :false-object :false)
               (json-error nil))))
    (unless (and (consp env) (consp (car env)))
      (genetics--error 'genetics-genome-error
                       "`%s' did not print a JSON envelope: %s"
                       (genetics-genome-argv-string argv)
                       (truncate-string-to-width (string-trim text) 200)))
    (unless (equal (alist-get 'schema env) genetics-genome-schema)
      (genetics--error 'genetics-genome-error
                       "`%s' answered schema %S, expected %S"
                       (genetics-genome-argv-string argv)
                       (alist-get 'schema env) genetics-genome-schema))
    env))

(defun genetics-genome-run (argv kind)
  "Run genome-cli ARGV synchronously and return its envelope alist.
The envelope must have kind KIND (a string).  A nonzero exit status or
an error envelope signals `genetics-genome-error' with genome-cli's code
and message."
  (let ((exe (executable-find (car argv))))
    (unless exe
      (genetics--error 'genetics-genome-missing
                       "genome-cli executable %S not found; install genome-cli or set `genetics-genome-executable' (or use `genetics-source-native')"
                       (car argv)))
    (let ((errfile (make-temp-file "genetics-genome-err")))
      (unwind-protect
          (with-temp-buffer
            (let* ((coding-system-for-read 'utf-8)
                   (status (apply #'call-process exe nil
                                  (list (current-buffer) errfile) nil
                                  (cdr argv)))
                   (out (buffer-string))
                   (err (with-temp-buffer
                          (insert-file-contents errfile)
                          (string-trim (buffer-string))))
                   (env (condition-case e
                            (genetics--genome-parse-envelope out argv)
                          (genetics-genome-error
                           (if (eq status 0)
                               (signal (car e) (cdr e))
                             (genetics--error
                              'genetics-genome-error
                              "`%s' failed (exit %s): %s"
                              (genetics-genome-argv-string argv) status
                              (if (string-empty-p err) "no output" err)))))))
              (when (or (not (eq status 0)) (eq (alist-get 'ok env) :false))
                (let ((e (alist-get 'error env)))
                  (genetics--error 'genetics-genome-error
                                   "genome %s failed (%s): %s"
                                   (cadr argv)
                                   (or (alist-get 'code e) (format "exit %s" status))
                                   (or (alist-get 'message e) err))))
              (unless (equal (alist-get 'kind env) kind)
                (genetics--error 'genetics-genome-error
                                 "`%s' answered kind %S, expected %S"
                                 (genetics-genome-argv-string argv)
                                 (alist-get 'kind env) kind))
              env))
        (delete-file errfile)))))

;;;; Mapping genome/v1 onto kits and records

(defun genetics--genome-symbol (value)
  "Return string VALUE as a symbol, or nil for null/empty.
Underscores become hyphens, so absent_means_ref and absent-means-ref
are the same."
  (and (stringp value) (not (string-empty-p value))
       (intern (string-replace "_" "-" value))))

(defun genetics--genome-format (source-format)
  "Return the kit format symbol for genome/v1 SOURCE-FORMAT."
  (genetics--genome-symbol source-format))

(defun genetics-genome-record->snp (rec)
  "Return a `genetics-snp' for genome/v1 genotype record REC.
Return nil when its call_source is \"missing\" (site not in the kit)."
  (let ((source (alist-get 'call_source rec)))
    (unless (equal source "missing")
      (let* ((chrom (genetics-normalize-chrom (format "%s" (alist-get 'chrom rec))))
             (pos (alist-get 'pos rec))
             (gt (alist-get 'genotype rec))
             (alt (alist-get 'alt rec)))
        (genetics-snp-create
         (or (alist-get 'rsid rec) (format "%s:%d" chrom pos))
         chrom pos
         (if (or (equal (alist-get 'zygosity rec) "no_call")
                 (not (stringp gt)) (string-empty-p gt))
             "--"
           (genetics--norm-genotype gt))
         (alist-get 'ref rec)
         (cond ((stringp alt) alt)
               (alt (mapconcat #'identity alt ",")))
         (when (equal source "inferred_ref") 'inferred-ref))))))

(defun genetics--genome-stats (summary)
  "Return the kit statistics plist for genome/v1 SUMMARY (one data item)."
  (let ((sex (alist-get 'sex summary))
        (by (alist-get 'by_chrom summary)))
    (list :total (alist-get 'records summary)
          :nocalls (alist-get 'no_calls summary)
          :het (alist-get 'het summary)
          :hom (+ (or (alist-get 'hom_alt summary) 0)
                  (or (alist-get 'hom_ref summary) 0))
          :hom-ref (alist-get 'hom_ref summary)
          :hemi (alist-get 'hemizygous summary)
          :chrom-counts
          (let ((pairs (mapcar (lambda (c)
                                 (cons (if (eq (car c) 'other_contigs)
                                           "other_contigs"
                                         (genetics-normalize-chrom
                                          (symbol-name (car c))))
                                       (cdr c)))
                               by)))
            (mapcar (lambda (c) (assoc c pairs))
                    (genetics--sort-chroms (mapcar #'car pairs))))
          :sex (when sex
                 (list :call (let ((c (alist-get 'call sex)))
                               (if (member c '("male" "female")) c "uncertain"))
                       :x-het-rate (alist-get 'x_het_rate sex)
                       :y-call-rate (alist-get 'y_call_rate sex)
                       :method (alist-get 'method sex)))
          :caveats (alist-get 'caveats summary))))

(defun genetics-genome-kit (kit-json summary-json file &optional warnings)
  "Build a `genetics-kit' from genome/v1 KIT-JSON and SUMMARY-JSON.
FILE is the file that was opened; WARNINGS are envelope warnings."
  (let* ((format (genetics--genome-format (alist-get 'source_format kit-json)))
         (build (genetics-build-key (alist-get 'build kit-json)))
         (note-format (if (memq format '(gvcf fastq-derived)) 'vcf format)))
    (genetics-kit--create
     :name (genetics--base-name file)
     :file file :format format :build build
     :chip (or (alist-get 'chip kit-json) "unknown")
     :sample (alist-get 'sample kit-json)
     :assay (genetics--genome-symbol (alist-get 'assay kit-json))
     :ref-calls (or (genetics--genome-symbol (alist-get 'ref_calls kit-json))
                    'unknown)
     :has-rsids (pcase (alist-get 'has_rsids kit-json)
                  (:false 'none) ('nil nil) (_ t))
     :strand-note (concat
                   (genetics--strand-note
                    note-format build
                    (equal (alist-get 'build_evidence kit-json) "assumed"))
                   " Parsed by genome-cli.")
     :stats (genetics--genome-stats summary-json)
     :backend #'genetics--genome-backend
     :backend-id (alist-get 'id kit-json)
     :warnings warnings
     :memo (make-hash-table :test 'equal))))

;;;; Backend

(defun genetics--genome-query (kit chrom offset &optional start end)
  "Return one page of records of KIT on CHROM from OFFSET.
START and END restrict positions.  The value is (SNPS . RAW-COUNT)."
  (let* ((env (genetics-genome-run
               (genetics-genome-argv 'query (genetics-kit-backend-id kit)
                                     chrom offset start end)
               "genotypes"))
         (data (alist-get 'data env)))
    (cons (delq nil (mapcar #'genetics-genome-record->snp data))
          (length data))))

(defun genetics--genome-walk (kit chrom fn &optional start end)
  "Call FN on each record of KIT on CHROM, page by page.
START and END restrict positions.  FN may exit non-locally."
  (let ((offset 0) (more t))
    (while more
      (pcase-let ((`(,snps . ,n)
                   (genetics--genome-query kit chrom offset start end)))
        (mapc fn snps)
        (setq offset (+ offset n)
              more (>= n genetics-genome-page-size))))))

(defun genetics--genome-get (kit rsid)
  "Return the call of RSID in genome-cli KIT (memoized), or nil."
  (let* ((memo (genetics-kit-memo kit))
         (hit (gethash rsid memo 'none)))
    (if (not (eq hit 'none))
        hit
      (let* ((env (genetics-genome-run
                   (genetics-genome-argv 'lookup (genetics-kit-backend-id kit)
                                         rsid)
                   "genotypes"))
             (snp (cl-some #'genetics-genome-record->snp
                           (alist-get 'data env))))
        (puthash rsid snp memo)))))

(defun genetics--genome-site-id (rec)
  "Return an id for genome/v1 genotype record or genotype REC."
  (cond ((stringp rec) rec)
        ((alist-get 'rsid rec))
        (t (format "%s:%s" (alist-get 'chrom rec) (alist-get 'pos rec)))))

(defun genetics--genome-pair (site)
  "Return (ID GENOTYPE-A GENOTYPE-B) for a genome/v1 discordant SITE.
SITE is an object with keys a and b (genotype records or strings), or a
list of two genotype records."
  (let* ((a (if (and (consp site) (consp (car site)) (symbolp (caar site)))
                (alist-get 'a site)
              (car site)))
         (b (if (and (consp site) (consp (car site)) (symbolp (caar site)))
                (alist-get 'b site)
              (cadr site)))
         (gt (lambda (r) (if (stringp r) r (alist-get 'genotype r)))))
    (list (genetics--genome-site-id
           (cond ((and (consp site) (consp (car site)) (symbolp (caar site))
                       (or (alist-get 'rsid site) (alist-get 'chrom site)))
                  site)
                 ((consp a) a)
                 ((consp b) b)
                 (t "?")))
          (funcall gt a) (funcall gt b))))

(defun genetics--genome-compare (a b)
  "Compare genome-cli kits A and B with `genome compare'."
  (let* ((env (genetics-genome-run
               (genetics-genome-argv 'compare (genetics-kit-backend-id a)
                                     (genetics-kit-backend-id b))
               "compare"))
         (r (car (alist-get 'data env)))
         (conc (alist-get 'concordance r))
         (disc (mapcar #'genetics--genome-pair (alist-get 'discordant_sites r))))
    (list :overlap (alist-get 'overlap r)
          :compared (+ (or (alist-get 'concordant r) 0)
                       (or (alist-get 'discordant r) 0))
          :concordant (alist-get 'concordant r)
          :discordant disc :complement nil
          :discordant-count (alist-get 'discordant r)
          :concordance (and conc (if (<= conc 1) (* 100.0 conc) conc))
          :build (alist-get 'build r)
          :warnings (alist-get 'warnings env)
          :build-mismatch nil)))

(defun genetics--genome-backend (op kit &rest args)
  "Answer kit operation OP for genome-cli KIT with ARGS.
OPs: `get' (RSID), `range' (CHROM START END), `map' (CHROM FN) and
`compare' (OTHER-KIT)."
  (pcase op
    ('get (genetics--genome-get kit (car args)))
    ('range (let (acc)
              (genetics--genome-walk kit (nth 0 args)
                                     (lambda (s) (push s acc))
                                     (nth 1 args) (nth 2 args))
              (nreverse acc)))
    ('map (genetics--genome-walk kit (car args) (cadr args)))
    ('compare (genetics--genome-compare kit (car args)))
    (_ (genetics--error 'genetics-genome-error "Unsupported operation %S" op))))

;;;; The genome-cli source

(defun genetics--genome-kits-for (file kits)
  "Entries of genome/v1 KITS envelope imported from FILE, newest first."
  (sort (seq-filter (lambda (k) (equal (alist-get 'source_path k) file))
                    (alist-get 'data kits))
        (lambda (a b) (string> (or (alist-get 'imported_at a) "")
                               (or (alist-get 'imported_at b) "")))))

(defun genetics--genome-fresh-p (kit-json file)
  "Non-nil when KIT-JSON was imported after FILE was last modified."
  (when-let* ((at (alist-get 'imported_at kit-json)))
    (not (time-less-p (date-to-time at)
                      (file-attribute-modification-time (file-attributes file))))))

(defun genetics-source-genome-cli (file &rest _args)
  "Open FILE through genome-cli and return a `genetics-kit'.
Reuses the kit genome-cli already holds for FILE when it was imported
after FILE last changed; otherwise runs `genome import FILE' (with
`--replace' when an older import exists), then `genome summary' (see
`genetics-source-genome-cli-explain').  Records stay in genome-cli and
are fetched on demand.  Keyword arguments of `genetics-parse-file' are
accepted and ignored."
  (let* ((file (expand-file-name file))
         (existing (genetics--genome-kits-for
                    file (genetics-genome-run (genetics-genome-argv 'kits) "kits")))
         (fresh (seq-find (lambda (k) (genetics--genome-fresh-p k file)) existing))
         (kits (if fresh
                   `((data ,fresh))
                 (genetics-genome-run (genetics-genome-argv 'import file (and existing t))
                                      "kits")))
         (kit-json (car (alist-get 'data kits))))
    (unless (alist-get 'id kit-json)
      (genetics--error 'genetics-genome-error
                       "genome import returned no kit for %s" file))
    (let* ((summary (genetics-genome-run
                     (genetics-genome-argv 'summary (alist-get 'id kit-json))
                     "summary"))
           (summary-json (car (alist-get 'data summary))))
      (genetics-genome-kit kit-json summary-json file
                           (append (alist-get 'warnings kits)
                                   (alist-get 'warnings summary))))))

(provide 'genetics-source)
;;; genetics-source.el ends here
