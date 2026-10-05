;;; genetics-org.el --- Org dynamic blocks for genetics.el -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Package-Requires: ((emacs "29.1"))
;; Keywords: data, tools, outlines
;; URL: https://github.com/davidawad/genetics.el

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; Three Org dynamic blocks, so a kit can be reported on from any Org
;; document (health-charts.el's report templates call them by name):
;;
;;   #+BEGIN: genetics-summary :file "~/dna/genome.txt"
;;   #+END:
;;   #+BEGIN: genetics-hits :kit "genome" :min-magnitude 2 :genes ("MTHFR")
;;   #+END:
;;   #+BEGIN: genetics-apoe :kit "genome"
;;   #+END:
;;
;; Update them with C-c C-x C-u (`org-dblock-update') or C-u C-c C-x C-u
;; for every block in the buffer.
;;
;; :kit names a loaded kit.  :file is a genotype file; a kit already
;; loaded from it is reused, otherwise it is opened through
;; `genetics-source-function' (genome-cli reuses its own import) and
;; registered without showing any buffer.  Relative :file names are
;; resolved from the Org file's directory.  With neither, the only loaded
;; kit is used.
;;
;; Every block ends with the "informational only, not medical advice"
;; line.  A failure never breaks the document: it is written as Org
;; comment lines with the error and what to do about it.  Each block has
;; a pure explain twin (`genetics-org-summary-explain' and friends, or
;; M-x `genetics-org-explain-block' on a block) that says what the block
;; would read and insert without opening or running anything.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'org)
(require 'org-table)
(require 'genetics-core)
(require 'genetics-source)
(require 'genetics-stats)
(require 'genetics-annotate)

(defconst genetics-org-disclaimer
  "/Informational only, not medical advice./ Consumer genotyping and sequencing are not diagnostic; confirm any finding with a clinical-grade test and a clinician."
  "Line ending every genetics dynamic block, including failed ones.")

(defconst genetics-org-blocks '("genetics-summary" "genetics-hits" "genetics-apoe")
  "Names of the dynamic blocks defined by genetics-org.el.")

;;;; Kit resolution

(defun genetics-org--kit-for-file (file)
  "Return the loaded kit read from FILE (expanded), or nil."
  (cl-find-if (lambda (k) (genetics--same-file-p (genetics-kit-file k) file))
              genetics-loaded-kits))

(defun genetics-org--param-file (params)
  "Return the expanded :file of PARAMS, or nil."
  (when-let* ((f (plist-get params :file)))
    (expand-file-name (format "%s" f))))

(defun genetics-org--param-kit (params)
  "Return the :kit name of PARAMS as a string, or nil."
  (when-let* ((k (plist-get params :kit)))
    (format "%s" k)))

(defun genetics-org-kit (params)
  "Return the kit designated by dynamic block PARAMS, opening it if needed.
:kit is the name of a loaded kit; :file a genotype file, reused when a
kit is loaded from it and otherwise opened with `genetics-source-open'
and registered (no buffer is shown).  With neither, the only loaded kit
is used.  Signal `genetics-no-kit' when no kit can be chosen."
  (let ((name (genetics-org--param-kit params))
        (file (genetics-org--param-file params)))
    (cond
     (name (genetics-find-kit name))
     (file
      (or (genetics-org--kit-for-file file)
          (progn
            (unless (file-readable-p file)
              (genetics--error 'genetics-file-error "Cannot read %s"
                               (abbreviate-file-name file)))
            (let ((kit (genetics-source-open file)))
              (setf (genetics-kit-name kit)
                    (genetics--unique-name (genetics-kit-name kit)
                                           (genetics-kit-file kit)))
              (genetics-register-kit kit)))))
     ((null genetics-loaded-kits)
      (genetics--error 'genetics-no-kit
                       "No kit loaded and the block has no :kit or :file"))
     ((null (cdr genetics-loaded-kits)) (car genetics-loaded-kits))
     (t (genetics--error 'genetics-no-kit
                         "%d kits are loaded; say which with :kit or :file"
                         (length genetics-loaded-kits))))))

(defun genetics-org-kit-explain (params)
  "Return a sentence saying which kit `genetics-org-kit' would use for PARAMS.
Pure: nothing is opened or run."
  (let ((name (genetics-org--param-kit params))
        (file (genetics-org--param-file params)))
    (cond
     (name
      (if (cl-find name genetics-loaded-kits :key #'genetics-kit-name :test #'equal)
          (format "Use the loaded kit %S." name)
        (format "Kit %S is not loaded; the block would fail (load it with M-x genetics-open)." name)))
     (file
      (let ((kit (genetics-org--kit-for-file file)))
        (cond
         (kit (format "Reuse the kit %S already loaded from %s." (genetics-kit-name kit)
                      (abbreviate-file-name file)))
         ((not (file-readable-p file))
          (format "%s is not readable; the block would fail." (abbreviate-file-name file)))
         ((or (eq genetics-source-function #'genetics-source-genome-cli)
              (and (eq genetics-source-function #'genetics-source-auto)
                   (genetics-genome-available-p)))
          (format "Open %s with genome-cli: `%s' (reusing genome-cli's import when it is newer than the file, otherwise `%s'), then `genome summary ID --format json'; register the kit without showing a buffer."
                  (abbreviate-file-name file)
                  (genetics-genome-argv-string (genetics-genome-argv 'kits))
                  (genetics-source-genome-cli-explain file)))
         ((eq genetics-source-function #'genetics-source-auto)
          (format "Parse %s with the Emacs Lisp parser (genome-cli is not installed) and register the kit without showing a buffer."
                  (abbreviate-file-name file)))
         (t (format "Open %s with `%s' and register the kit without showing a buffer."
                    (abbreviate-file-name file) genetics-source-function)))))
     ((null genetics-loaded-kits)
      "No :kit or :file and no kit is loaded; the block would fail.")
     ((null (cdr genetics-loaded-kits))
      (format "Use the only loaded kit, %S." (genetics-kit-name (car genetics-loaded-kits))))
     (t "Several kits are loaded and the block names none; it would fail."))))

;;;; Failures

(defun genetics-org-runbook (err)
  "Return what the user should do about ERR, a condition object."
  (pcase (car err)
    ('genetics-no-kit
     "Load the kit with M-x genetics-open, or give the block :file \"path/to/kit\" (or :kit with a loaded kit's name).")
    ((or 'genetics-file-error 'file-error 'file-missing)
     "Check that :file names an existing, readable genotype file; relative names are resolved from this Org file's directory.")
    ((or 'genetics-unknown-format 'genetics-parse-error)
     "Readable formats: 23andMe, AncestryDNA, MyHeritage/FTDNA CSV, VCF and VCF.gz (see README, \"Supported formats\").")
    ((or 'genetics-unsupported-file 'genetics-fastq-file)
     "FASTQ/BAM/CRAM hold reads, not genotypes: make a VCF with M-x genetics-fastq-plan (genome-cli), then point :file at the VCF.")
    ('genetics-gzip-error "Install gzip (or use an Emacs built with zlib), or decompress the VCF first.")
    ('genetics-genome-missing
     "Install genome-cli, set `genetics-genome-executable', or set `genetics-source-function' to `genetics-source-native'.")
    ('genetics-genome-error
     "Run the command shown by M-x genetics-source-genome-cli-explain in a shell to see genome-cli's own error.")
    ('genetics-annotation-error
     "Fix the annotation file named above (see README, \"Annotation files\"), then M-x genetics-reload-annotations.")
    (_ "M-x genetics-org-explain-block on this block shows what it tries to read.")))

(defun genetics-org-error-text (block err)
  "Return Org comment lines reporting ERR in dynamic block BLOCK.
The disclaimer line follows, so even a failed block carries it."
  (let ((msg (if (and (get (car err) 'error-conditions)
                       (memq 'genetics-error (get (car err) 'error-conditions))
                       (stringp (cadr err)))
                  (format "%s: %s" (get (car err) 'error-message) (cadr err))
                (error-message-string err))))
    (concat (format "# %s failed: %s\n" block
                    (replace-regexp-in-string "\n" " " msg))
            (format "# What to do: %s\n\n" (genetics-org-runbook err))
            genetics-org-disclaimer)))

(defun genetics-org-block-string (block params fn)
  "Return the text of dynamic BLOCK for PARAMS rendered by FN, never signalling.
FN takes a kit and PARAMS and returns Org text.  Any error becomes the
comment lines of `genetics-org-error-text'."
  (condition-case err
      (concat (string-trim-right (funcall fn (genetics-org-kit params) params))
              "\n\n" genetics-org-disclaimer)
    (error (genetics-org-error-text block err))))

;;;; Shared formatting

(defun genetics-org--cell (value)
  "Return VALUE formatted as a safe Org table cell."
  (let ((s (if value (format "%s" value) "")))
    (string-trim (replace-regexp-in-string
                  "[ \t\n]+" " " (replace-regexp-in-string "|" "/" s)))))

(defun genetics-org--align (text)
  "Return Org TEXT with every table aligned."
  (with-temp-buffer
    (insert text)
    (delay-mode-hooks (org-mode))
    ;; Fontify so hidden link targets do not count towards column widths.
    (font-lock-ensure)
    (goto-char (point-min))
    (while (re-search-forward org-table-line-regexp nil t)
      (org-table-align)
      (goto-char (org-table-end)))
    (buffer-string)))

(defun genetics-org--first-sentence (text)
  "Return the first sentence of TEXT (all of it if there is no full stop)."
  (if (and text (string-match "\\`\\(.+?[.!?]\\)\\(?: \\|\\'\\)" text))
      (match-string 1 text)
    text))

;;;; genetics-summary

(defun genetics-org-summary-string (kit &optional _params)
  "Return the Org section body summarising KIT.
Pure given KIT: format, source, assay, build, records, no-call rate,
inferred sex and caveats."
  (let* ((s (genetics-kit-stats kit))
         (nocalls (plist-get s :nocalls))
         (rate (genetics-nocall-rate kit)))
    (with-temp-buffer
      (insert (format "- Kit :: %s\n" (genetics-kit-name kit)))
      (insert (format "- File :: =%s=\n" (abbreviate-file-name (genetics-kit-file kit))))
      (insert (format "- Format :: %s\n" (genetics--format-label (genetics-kit-format kit))))
      (insert (format "- Source :: %s\n"
                      (if (genetics-kit-backend kit)
                          (format "genome-cli (kit %s)" (genetics-kit-backend-id kit))
                        "Emacs Lisp parser")))
      (insert (format "- Assay :: %s%s\n" (or (genetics-kit-assay kit) "unknown")
                      (pcase (genetics-kit-ref-calls kit)
                        ('absent-means-ref ", variant sites only (absent = reference, inferred)")
                        ('explicit ", every assayed site listed")
                        (_ ""))))
      (insert (format "- Build :: GRCh%s\n" (or (genetics-kit-build kit) "unknown")))
      (when (and (genetics-kit-chip kit) (not (equal (genetics-kit-chip kit) "unknown")))
        (insert (format "- Chip :: %s\n" (genetics-kit-chip kit))))
      (when (genetics-kit-sample kit)
        (insert (format "- Sample :: %s\n" (genetics-kit-sample kit))))
      (insert (format "- Records :: %d\n" (or (plist-get s :total) 0)))
      (insert (format "- No-call rate :: %s\n"
                      (if (and nocalls rate)
                          (format "%.2f%% (%d no-calls)" (* 100 rate) nocalls)
                        "not computed (records are read on demand)")))
      (insert (format "- Inferred sex :: %s\n"
                      (if (and (genetics-kit-lazy kit) (not (plist-get s :sex)))
                          "not computed (offset-indexed)"
                        (genetics-sex-description kit))))
      (insert "- Caveats ::\n")
      (dolist (c (genetics-kit-caveats kit))
        (insert (format "  - %s\n" c)))
      (buffer-string))))

(defun genetics-org-summary-explain (params)
  "Return what the genetics-summary block would do for PARAMS.
Pure twin of `org-dblock-write:genetics-summary': nothing is opened or run."
  (concat "genetics-summary: " (genetics-org-kit-explain params)
          " Insert format, source, assay, build, records, no-call rate, inferred sex and caveats, then the not-medical-advice line."))

;;;###autoload
(defun org-dblock-write:genetics-summary (params)
  "Write the genetics-summary dynamic block for PARAMS (:kit or :file).
See `genetics-org-summary-string'; failures become Org comment lines."
  (insert (genetics-org-block-string "genetics-summary" params
                                     #'genetics-org-summary-string)))

;;;; genetics-hits

(defun genetics-org--genes (genes)
  "Normalise the :genes parameter GENES to a list of upper-case strings.
GENES may be a list of strings or symbols, or one string of names
separated by spaces or commas."
  (mapcar #'upcase
          (cond ((null genes) nil)
                ((listp genes) (mapcar (lambda (g) (format "%s" g)) genes))
                (t (split-string (format "%s" genes) "[ ,]+" t)))))

(defun genetics-org--number (value)
  "Return VALUE (a number, or a string or symbol naming one) as a number.
Return nil for nil or anything that is not a number."
  (cond ((numberp value) value)
        ((and value (string-match-p "\\`[0-9]+\\(\\.[0-9]+\\)?\\'" (format "%s" value)))
         (string-to-number (format "%s" value)))))

(defun genetics-org--magnitude (ann)
  "Return the magnitude of ANN as a number, or nil."
  (genetics-org--number (genetics-annotation-magnitude ann)))

(defun genetics-org-filter-hits (hits min-magnitude genes)
  "Return the (ANNOTATION . SNP) HITS passing MIN-MAGNITUDE and GENES.
A hit without a magnitude fails any MIN-MAGNITUDE."
  (let ((genes (genetics-org--genes genes)))
    (cl-remove-if-not
     (lambda (h)
       (let ((ann (car h)))
         (and (or (null min-magnitude)
                  (let ((m (genetics-org--magnitude ann)))
                    (and m (>= m min-magnitude))))
              (or (null genes)
                  (member (upcase (or (genetics-annotation-gene ann) "")) genes)))))
     hits)))

(defun genetics-org--hit-row (ann snp)
  "Return the table cells for annotation ANN and call SNP."
  (let* ((a (genetics-assess-snp ann snp))
         (copies (plist-get a :copies))
         (risk (genetics-annotation-risk-allele ann))
         (url (genetics-annotation-url ann)))
    (list (genetics-annotation-rsid ann)
          (genetics-annotation-gene ann)
          (genetics-snp-genotype snp)
          (cond ((null copies) "n/a")
                (risk (format "%s (%s)" copies risk))
                (t copies))
          (if (genetics-snp-inferred-p snp) "inferred ref" "observed")
          (genetics-annotation-magnitude ann)
          (genetics-org--first-sentence (genetics-annotation-effect ann))
          (if url (format "[[%s][source]]" url) ""))))

(defun genetics-org-hits-string (kit &optional params)
  "Return an Org table of the curated annotation hits in KIT.
PARAMS may hold :min-magnitude (a number), :genes (see
`genetics-org--genes') and :effect-width (a number: add an Org width
cookie so the effect column can be shrunk, see `org-table-shrink').
Strand flags and inferred calls are explained under the table."
  (let* ((min (genetics-org--number (plist-get params :min-magnitude)))
         (width (genetics-org--number (plist-get params :effect-width)))
         (genes (plist-get params :genes))
         (all (genetics-annotated-snps kit))
         (hits (genetics-org-filter-hits all min genes))
         notes)
    (with-temp-buffer
      (if (null hits)
          (insert (if all
                      (format "No annotated SNP in this kit passes the filter (%d before filtering).\n"
                              (length all))
                    "No annotated SNPs were found in this kit.\n"))
        (insert "| rsid | gene | genotype | risk-allele copies | call source | magnitude | effect | source |\n|-\n")
        (when (natnump width)
          (insert (format "| | | | | | | <%d> | |\n" width)))
        (dolist (h hits)
          (insert "| " (mapconcat #'genetics-org--cell
                                  (genetics-org--hit-row (car h) (cdr h)) " | ")
                  " |\n")
          (let* ((a (genetics-assess-snp (car h) (cdr h)))
                 (flag (genetics-flag-text (plist-get a :flag)
                                           (plist-get a :flipped-copies))))
            (when (and flag (not (eq (plist-get a :flag) 'no-call)))
              (push (format "- %s :: %s" (genetics-annotation-rsid (car h)) flag)
                    notes))))
        (when (cl-some (lambda (h) (genetics-snp-inferred-p (cdr h))) hits)
          (insert "\nCall source \"inferred ref\": " genetics-inferred-ref-text "\n"))
        (when notes
          (insert "\nStrand notes:\n" (string-join (nreverse notes) "\n") "\n"))
        (insert "\nRisk alleles are on the + strand; a genotype is never flipped automatically.\n"))
      (when (or min genes)
        (insert (format "\nFilter: %s.\n"
                        (string-join
                         (delq nil (list (and min (format "magnitude >= %s" min))
                                         (and genes (format "genes %s"
                                                            (string-join (genetics-org--genes genes) ", ")))))
                         "; "))))
      (genetics-org--align (buffer-string)))))

(defun genetics-org-hits-explain (params)
  "Return what the genetics-hits block would do for PARAMS.
Pure twin of `org-dblock-write:genetics-hits': nothing is opened or run."
  (let ((min (genetics-org--number (plist-get params :min-magnitude)))
        (genes (genetics-org--genes (plist-get params :genes))))
    (format "genetics-hits: %s Tabulate the curated annotation hits%s%s (rsid, gene, genotype, risk-allele copies, call source, magnitude, effect, source link) from %s, then the not-medical-advice line."
            (genetics-org-kit-explain params)
            (if min (format " with magnitude >= %s" min) "")
            (if genes (format " in %s" (string-join genes ", ")) "")
            (string-join (mapcar #'abbreviate-file-name genetics-annotation-files) ", "))))

;;;###autoload
(defun org-dblock-write:genetics-hits (params)
  "Write the genetics-hits dynamic block for PARAMS.
PARAMS: :kit or :file, :min-magnitude, :genes, :effect-width.  See
`genetics-org-hits-string'; failures become Org comment lines."
  (insert (genetics-org-block-string "genetics-hits" params
                                     #'genetics-org-hits-string)))

;;;; genetics-apoe

(defun genetics-org-apoe-string (kit &optional _params)
  "Return an Org paragraph interpreting the APOE haplotype of KIT.
Carries the report's caveats: inferred calls are labelled, the strand
convention is stated and phase ambiguity is explained."
  (let ((apoe (genetics-apoe-for-kit kit)))
    (if (null apoe)
        "APOE could not be assessed: rs429358 and rs7412 were not both present in this kit.\n"
      (let ((calls (mapcar (lambda (rsid)
                             (let ((snp (genetics-kit-resolve kit rsid)))
                               (format "%s %s" rsid (genetics-snp-genotype-label snp))))
                           genetics-apoe-snps)))
        (concat
         (format "*APOE: %s* (status: %s). " (or (plist-get apoe :diplotype) "undetermined")
                 (plist-get apoe :status))
         (plist-get apoe :description)
         (format " Calls: %s." (string-join calls ", "))
         " Method: rs429358 T/C and rs7412 C/T read on the + strand (e2 = T+T, e3 = T+C, e4 = C+C); genotypes are never flipped automatically, and an A/G call is reported as a possible strand flip."
         (when (plist-get apoe :inferred)
           (concat " " genetics-inferred-ref-text))
         " APOE e4 is a risk factor, not a diagnosis; most carriers never develop Alzheimer's disease.\n")))))

(defun genetics-org-apoe-explain (params)
  "Return what the genetics-apoe block would do for PARAMS.
Pure twin of `org-dblock-write:genetics-apoe': nothing is opened or run."
  (concat "genetics-apoe: " (genetics-org-kit-explain params)
          " Resolve rs429358 and rs7412 (by position on the kit's build when it has no rsids; inferred reference calls are labelled), write the APOE diplotype paragraph with strand notes, then the not-medical-advice line."))

;;;###autoload
(defun org-dblock-write:genetics-apoe (params)
  "Write the genetics-apoe dynamic block for PARAMS (:kit or :file).
See `genetics-org-apoe-string'; failures become Org comment lines."
  (insert (genetics-org-block-string "genetics-apoe" params
                                     #'genetics-org-apoe-string)))

;;;; Commands

(defun genetics-org-explain (block params)
  "Return the explain text of dynamic BLOCK (a name) for PARAMS."
  (pcase block
    ("genetics-summary" (genetics-org-summary-explain params))
    ("genetics-hits" (genetics-org-hits-explain params))
    ("genetics-apoe" (genetics-org-apoe-explain params))
    (_ (genetics--error 'genetics-error "%s is not a genetics dynamic block" block))))

;;;###autoload
(defun genetics-org-explain-block ()
  "Show what the genetics dynamic block at point would read and insert.
Nothing is opened or run.  Returns the explanation."
  (interactive)
  (save-excursion
    (unless (re-search-backward org-dblock-start-re nil t)
      (genetics--error 'genetics-error "Not inside a dynamic block"))
    (let* ((params (genetics-org--block-params))
           (text (genetics-org-explain (plist-get params :name) params)))
      (when (called-interactively-p 'interactive)
        (message "%s" text))
      text)))

(defun genetics-org--block-params ()
  "Return the parameter plist of the dynamic block starting at point."
  (save-excursion
    (beginning-of-line)
    (unless (looking-at org-dblock-start-re)
      (genetics--error 'genetics-error "Not at a dynamic block"))
    (let ((name (org-no-properties (match-string 1)))
          (args (match-string 3)))
      (append (list :name name)
              (and args (read (concat "(" args ")")))))))

;;;###autoload
(defun genetics-org-insert-block (block &optional kit)
  "Insert an empty dynamic BLOCK for KIT (a kit name) at point and fill it."
  (interactive
   (list (completing-read "Block: " genetics-org-blocks nil t)
         (when genetics-loaded-kits
           (completing-read "Kit: " (mapcar #'genetics-kit-name genetics-loaded-kits)
                            nil t))))
  (unless (bolp) (insert "\n"))
  (insert (format "#+BEGIN: %s%s\n#+END:\n" block
                  (if kit (format " :kit %S" kit) "")))
  (forward-line -2)
  (org-update-dblock))

(provide 'genetics-org)
;;; genetics-org.el ends here
