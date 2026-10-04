;;; genetics-core.el --- Shared definitions for genetics.el -*- lexical-binding: t; -*-

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

;; Customization, error types, data structures and small helpers shared by
;; every genetics.el module.  Nothing here touches the network.

;;; Code:

(require 'cl-lib)
(require 'subr-x)

(defgroup genetics nil
  "Read and explore consumer genetics raw-data files."
  :group 'tools
  :prefix "genetics-")

(defconst genetics--package-dir
  (file-name-directory (or load-file-name buffer-file-name default-directory))
  "Directory the genetics package was loaded from.")

;;;; Customization

(defcustom genetics-data-directory
  "~/Documents/Genetics/"
  "Default directory offered when prompting for a raw-data file."
  :type 'directory)

(defcustom genetics-cache-directory (locate-user-emacs-file "genetics-cache/")
  "Directory for parse caches, decompressed VCFs and SNPedia answers."
  :type 'directory)

(defcustom genetics-use-cache t
  "Non-nil means cache parsed kits on disk as `.eld' files."
  :type 'boolean)

(defcustom genetics-vcf-eager-limit (* 200 1024 1024)
  "VCF files larger than this many bytes are offset-indexed, not loaded."
  :type 'integer)

(defcustom genetics-chunk-size (* 4 1024 1024)
  "Bytes read at a time when streaming a large file."
  :type 'integer)

(defcustom genetics-browse-limit 5000
  "Maximum number of rows shown in the browser buffer."
  :type 'integer)

(defcustom genetics-annotation-files
  (list (expand-file-name "annotations/genetics-curated.json"
                          genetics--package-dir))
  "List of annotation files (JSON or Org); later files override earlier."
  :type '(repeat file))

(defcustom genetics-vcf-ref-calls 'auto
  "How absent sites of a VCF parsed in Emacs are treated.
`auto' treats a VCF as variant-only whole-genome data (`absent-means-ref')
when it has no explicit homozygous-reference calls and at least
`genetics-wgs-min-records' records, and as `explicit' when it is a gVCF;
otherwise `unknown'.  `absent-means-ref' and `unknown' force that
answer.  Only `absent-means-ref' lets curated sites that are missing from
the file be shown as inferred homozygous reference."
  :type '(choice (const :tag "Detect" auto)
                 (const :tag "Variant-only WGS (absent = reference)"
                        absent-means-ref)
                 (const :tag "Never infer" unknown)))

(defcustom genetics-wgs-min-records 1000000
  "Minimum VCF record count for `auto' to treat a VCF as whole-genome."
  :type 'integer)

(defcustom genetics-snpedia-enabled nil
  "Non-nil allows SNPedia lookups, which send the rsid over the network."
  :type 'boolean)

;;;; Errors

(define-error 'genetics-error "Genetics error")
(define-error 'genetics-file-error "Genetics file error" 'genetics-error)
(define-error 'genetics-unknown-format "Unrecognized genetics file format"
  'genetics-error)
(define-error 'genetics-parse-error "Cannot parse genetics data"
  'genetics-error)
(define-error 'genetics-gzip-error "Cannot decompress with gzip"
  'genetics-error)
(define-error 'genetics-unsupported-file
  "Raw reads cannot be opened as genotypes" 'genetics-error)
(define-error 'genetics-fastq-file "FASTQ reads cannot be opened directly"
  'genetics-unsupported-file)
(define-error 'genetics-genome-error "genome-cli failed" 'genetics-error)
(define-error 'genetics-genome-missing "genome-cli executable not found"
  'genetics-genome-error)
(define-error 'genetics-no-kit "No genetics kit available" 'genetics-error)
(define-error 'genetics-annotation-error "Bad genetics annotation"
  'genetics-error)
(define-error 'genetics-snpedia-error "SNPedia lookup failed" 'genetics-error)
(define-error 'genetics-snpedia-disabled
  "SNPedia lookups are disabled (set `genetics-snpedia-enabled')"
  'genetics-snpedia-error)
(define-error 'genetics-snpedia-declined "SNPedia lookup not confirmed"
  'genetics-snpedia-error)

(defun genetics--error (type fmt &rest args)
  "Signal error TYPE with a message built from FMT and ARGS."
  (signal type (list (apply #'format fmt args))))

;;;; Records

(cl-defstruct (genetics-snp (:constructor genetics-snp-create
                                          (rsid chrom pos genotype
                                                &optional ref alt call-source))
                            (:copier nil))
  "One genotype call.  GENOTYPE is upper case letters or `--' for no-call.
CALL-SOURCE is nil (observed in the file) or `inferred-ref' (absent from
a variant-only whole-genome VCF and assumed homozygous reference)."
  rsid chrom pos genotype ref alt call-source)

(cl-defstruct (genetics-kit (:constructor genetics-kit--create)
                            (:copier nil))
  "A loaded raw-data file.
ASSAY is `array', `wgs', `wes', `panel' or `unknown'.  REF-CALLS says how
reference calls are represented: `explicit' (every assayed site is listed),
`absent-means-ref' (variant-only WGS: a covered site that is absent is
homozygous reference) or `unknown'.  BACKEND is nil for kits parsed in
Emacs, or a function implementing the genome-cli source (see
genetics-source.el); BACKEND-ID is the kit id in that source.
HAS-RSIDS is t, `none' or nil (unknown)."
  name file format build chip strand-note sample
  table chroms stats lazy index ranges data-file line-parser blocks
  assay ref-calls has-rsids backend backend-id warnings memo)

(defun genetics-snp-inferred-p (snp)
  "Return non-nil if SNP was inferred rather than observed."
  (eq (genetics-snp-call-source snp) 'inferred-ref))

(defun genetics-snp-genotype-label (snp)
  "Return the display genotype of SNP; an inferred call is marked."
  (if (genetics-snp-inferred-p snp)
      (format "%s (inferred ref)" (genetics-snp-genotype snp))
    (genetics-snp-genotype snp)))

(defconst genetics-inferred-ref-text
  "Inferred, not observed: this site is absent from a variant-only whole-genome VCF, so it is assumed homozygous for the reference allele of the kit's build. A site that was not covered by sequencing would look the same."
  "Explanation attached to every inferred homozygous-reference call.")

(defvar genetics-loaded-kits nil
  "List of loaded `genetics-kit' objects, most recent first.")

(defvar-local genetics--buffer-kit nil
  "Kit shown by the current genetics buffer.")

;;;; Chromosomes

(defconst genetics-chromosome-order
  '("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16"
    "17" "18" "19" "20" "21" "22" "X" "Y" "XY" "MT")
  "Canonical chromosome order.")

(defvar genetics--chrom-table (make-hash-table :test 'equal)
  "Intern table so equal chromosome strings are `eq'.")

(defun genetics--intern-chrom (chrom)
  "Return the shared copy of string CHROM."
  (or (gethash chrom genetics--chrom-table)
      (puthash chrom chrom genetics--chrom-table)))

(defun genetics-normalize-chrom (chrom &optional numeric-sex)
  "Normalize CHROM (drop chr prefix, M to MT).
With NUMERIC-SEX map 23, 24, 25, 26 to X, Y, XY, MT (AncestryDNA)."
  (let ((c (upcase (string-trim chrom))))
    (when (string-prefix-p "CHR" c)
      (setq c (substring c 3)))
    (genetics--intern-chrom
     (cond ((member c '("M" "MT" "MITO")) "MT")
           ((not numeric-sex) c)
           ((equal c "23") "X")
           ((equal c "24") "Y")
           ((equal c "25") "XY")
           ((equal c "26") "MT")
           (t c)))))

(defun genetics--chrom-rank (chrom)
  "Return a sort rank for CHROM."
  (or (cl-position chrom genetics-chromosome-order :test #'equal)
      (+ 100 (sxhash-equal chrom))))

(defun genetics-primary-chrom-p (chrom)
  "Return non-nil if CHROM is a primary assembly chromosome (1-22, X, Y, XY, MT)."
  (member chrom genetics-chromosome-order))

(defun genetics-fold-chrom-counts (counts)
  "Fold non-primary contigs of COUNTS, an alist (CHROM . N), into one row.
Primary chromosomes keep their rows in order; every alt, decoy, HLA or
unplaced contig is summed into a final row (\"other contigs (K)\" . N)
where K is the number of such contigs.  An existing \"other_contigs\"
entry (genome-cli summaries) is added to that row."
  (let ((primary nil) (k 0) (n 0) (seen nil))
    (dolist (c counts)
      (cond ((genetics-primary-chrom-p (car c)) (push c primary))
            ((equal (car c) "other_contigs")
             (setq seen t n (+ n (cdr c))))
            (t (setq k (1+ k) n (+ n (cdr c))))))
    (nconc (nreverse primary)
           (cond ((> k 0) (list (cons (format "other contigs (%d)" k) n)))
                 (seen (list (cons "other contigs" n)))))))

(defun genetics--sort-chroms (chroms)
  "Return CHROMS (list of strings) in canonical order."
  (sort (copy-sequence chroms)
        (lambda (a b)
          (let ((ra (genetics--chrom-rank a)) (rb (genetics--chrom-rank b)))
            (if (= ra rb) (string< a b) (< ra rb))))))

;;;; Genotypes

(defun genetics-no-call-p (genotype)
  "Return non-nil if GENOTYPE is a no-call."
  (or (null genotype) (string-empty-p genotype) (string= genotype "--")))

(defun genetics-alleles (genotype)
  "Return the list of allele strings in GENOTYPE (nil for a no-call)."
  (cond ((genetics-no-call-p genotype) nil)
        ((string-search "/" genotype) (split-string genotype "/" t))
        (t (mapcar #'char-to-string genotype))))

(defun genetics-genotype-key (genotype)
  "Return an order-insensitive key for GENOTYPE, or nil for a no-call.
A single (hemizygous) allele is treated as its homozygous form."
  (let ((a (genetics-alleles genotype)))
    (when a
      (when (null (cdr a))
        (setq a (list (car a) (car a))))
      (mapconcat #'identity (sort (copy-sequence a) #'string<) "/"))))

(defun genetics-zygosity (genotype)
  "Classify GENOTYPE as no-call, heterozygous, homozygous or hemizygous."
  (let ((a (genetics-alleles genotype)))
    (pcase (length a)
      (0 'no-call)
      (1 'hemizygous)
      (2 (if (string= (car a) (cadr a)) 'homozygous 'heterozygous))
      (_ 'other))))

(defun genetics-complement (allele)
  "Return the complement of single-base ALLELE, or nil if not A/C/G/T."
  (pcase allele
    ("A" "T") ("T" "A") ("C" "G") ("G" "C") (_ nil)))

;;;; Streaming large files

(defun genetics--stream-lines (file start end fn)
  "Call FN with each line of FILE between byte offsets START and END.
FN receives the line (without newline) and its starting byte offset.
END nil means end of file.  Reads `genetics-chunk-size' bytes at a time.
If FN returns the symbol `stop', streaming ends."
  (let* ((end (or end (file-attribute-size (file-attributes file))))
         (pos start)
         (line-start start)
         (partial nil)
         (done nil))
    (with-temp-buffer
      (set-buffer-multibyte nil)
      (while (and (not done) (< pos end))
        (erase-buffer)
        (insert-file-contents-literally
         file nil pos (min end (+ pos (max 1 genetics-chunk-size))))
        (let ((n (buffer-size)) (from (point-min)))
          (when (zerop n) (setq done t))
          (goto-char (point-min))
          (while (and (not done) (search-forward "\n" nil t))
            (let ((line (buffer-substring-no-properties from (1- (point)))))
              (when partial
                (setq line (concat partial line) partial nil))
              (when (eq (funcall fn line line-start) 'stop)
                (setq done t))
              (setq line-start (+ pos (- (point) (point-min))))
              (setq from (point))))
          (unless done
            (when (< from (point-max))
              (setq partial (concat partial (buffer-substring-no-properties
                                             from (point-max)))))
            (setq pos (+ pos n)))))
      (when (and partial (not done) (> (length partial) 0))
        (funcall fn partial line-start)))))

(defun genetics--read-line-at (file offset)
  "Return the line of FILE beginning at byte OFFSET."
  (let ((chunk 4096) (acc "") (pos offset) (result nil))
    (with-temp-buffer
      (set-buffer-multibyte nil)
      (while (not result)
        (erase-buffer)
        (insert-file-contents-literally file nil pos (+ pos chunk))
        (if (zerop (buffer-size))
            (setq result acc)
          (goto-char (point-min))
          (if (search-forward "\n" nil t)
              (setq result (concat acc (buffer-substring-no-properties
                                        (point-min) (1- (point)))))
            (setq acc (concat acc (buffer-string)) pos (+ pos (buffer-size)))))))
    result))

;;;; Kit access (works for eager and offset-indexed kits)

(defun genetics-kit-get (kit rsid)
  "Return the `genetics-snp' for RSID in KIT, or nil."
  (cond ((genetics-kit-backend kit)
         (funcall (genetics-kit-backend kit) 'get kit rsid))
        ((genetics-kit-lazy kit)
         (let ((off (gethash rsid (genetics-kit-index kit))))
           (when off
             (funcall (genetics-kit-line-parser kit)
                      (genetics--read-line-at (genetics-kit-data-file kit)
                                              off)))))
        (t (gethash rsid (genetics-kit-table kit)))))

(defun genetics-kit-ids (kit)
  "Return a hash table whose keys are the ids available in KIT.
Kits served by genome-cli return an empty table (ids stay in genome-cli)."
  (cond ((genetics-kit-backend kit) (make-hash-table :test 'equal))
        ((genetics-kit-lazy kit) (genetics-kit-index kit))
        (t (genetics-kit-table kit))))

(defun genetics-kit-snp-count (kit)
  "Return the number of records in KIT."
  (if (genetics-kit-backend kit)
      (plist-get (genetics-kit-stats kit) :total)
    (hash-table-count (genetics-kit-ids kit))))

(defun genetics-kit-chromosomes (kit)
  "Return the chromosomes present in KIT in canonical order."
  (cond ((genetics-kit-backend kit)
         (cl-remove-if-not #'genetics-primary-chrom-p
                           (mapcar #'car (plist-get (genetics-kit-stats kit)
                                                    :chrom-counts))))
        ((genetics-kit-lazy kit)
         (genetics--sort-chroms (mapcar #'car (genetics-kit-ranges kit))))
        (t (mapcar #'car (genetics-kit-chroms kit)))))

(defun genetics-kit-map-snps (kit fn &optional chrom)
  "Call FN on each SNP of KIT in chromosome/position order.
Only CHROM is visited when given.  FN returning `stop' ends the walk."
  (let ((chroms (if chrom (list chrom) (genetics-kit-chromosomes kit))))
    (catch 'genetics--stop
      (dolist (c chroms)
        (cond
         ((genetics-kit-backend kit)
          (funcall (genetics-kit-backend kit) 'map kit c
                   (lambda (snp)
                     (when (eq (funcall fn snp) 'stop)
                       (throw 'genetics--stop nil)))))
         ((genetics-kit-lazy kit)
          (genetics--lazy-map-chrom kit c fn))
         (t
          (let ((vec (cdr (assoc c (genetics-kit-chroms kit)))))
            (cl-loop for snp across vec
                     when (eq (funcall fn snp) 'stop)
                     do (throw 'genetics--stop nil)))))))))

(defun genetics--lazy-map-chrom (kit chrom fn)
  "Stream the records of CHROM in lazy KIT through FN."
  (let ((range (cdr (assoc chrom (genetics-kit-ranges kit))))
        (parser (genetics-kit-line-parser kit)))
    (when range
      (genetics--stream-lines
       (genetics-kit-data-file kit) (car range) (cdr range)
       (lambda (line _off)
         (unless (or (string-empty-p line) (eq (aref line 0) ?#))
           (let ((snp (funcall parser line)))
             (when (and snp (equal (genetics-snp-chrom snp) chrom)
                        (eq (funcall fn snp) 'stop))
               (throw 'genetics--stop nil)))))))))

;;;; Access by position

(defun genetics--lower-bound (vec pos key)
  "Return the first index of sorted VEC whose KEY value is >= POS."
  (let ((lo 0) (hi (length vec)))
    (while (< lo hi)
      (let ((mid (/ (+ lo hi) 2)))
        (if (< (funcall key (aref vec mid)) pos)
            (setq lo (1+ mid))
          (setq hi mid))))
    lo))

(defun genetics-kit-records-in (kit chrom start end)
  "Return the records of KIT on CHROM with START <= position <= END.
Eager kits use binary search; offset-indexed kits seek with their block
index and stream only the lines in between."
  (cond
   ((genetics-kit-backend kit)
    (funcall (genetics-kit-backend kit) 'range kit chrom start end))
   ((genetics-kit-lazy kit)
    (let* ((range (cdr (assoc chrom (genetics-kit-ranges kit))))
           (blocks (cdr (assoc chrom (genetics-kit-blocks kit))))
           (parser (genetics-kit-line-parser kit))
           (from (when range
                   (if (and blocks (> (length blocks) 0))
                       (let ((i (genetics--lower-bound blocks start #'car)))
                         ;; start one block early: equal positions may span
                         (if (> i 0) (cdr (aref blocks (1- i))) (car range)))
                     (car range))))
           (acc nil))
      (when range
        (genetics--stream-lines
         (genetics-kit-data-file kit) from (cdr range)
         (lambda (line _off)
           (unless (or (string-empty-p line) (eq (aref line 0) ?#))
             (let ((snp (funcall parser line)))
               (when (equal (genetics-snp-chrom snp) chrom)
                 (cond ((> (genetics-snp-pos snp) end) 'stop)
                       ((>= (genetics-snp-pos snp) start) (push snp acc)
                        nil))))))))
      (nreverse acc)))
   (t
    (let ((vec (cdr (assoc chrom (genetics-kit-chroms kit)))) (acc nil))
      (when vec
        (cl-loop for i from (genetics--lower-bound vec start #'genetics-snp-pos)
                 below (length vec)
                 for snp = (aref vec i)
                 while (<= (genetics-snp-pos snp) end)
                 do (push snp acc)))
      (nreverse acc)))))

(defun genetics--snv-record-p (snp)
  "Return non-nil if SNP is a single-base site (or has no REF)."
  (let ((ref (genetics-snp-ref snp)))
    (or (null ref) (= (length ref) 1))))

(defun genetics-kit-at (kit chrom pos &optional ref)
  "Return the observed single-base record of KIT at CHROM:POS, or nil.
Indels anchored at the position are ignored: their anchor base is the
reference, so they say nothing about a SNP there.  With REF, a record
whose REF matches is preferred."
  (let ((hits (cl-remove-if-not #'genetics--snv-record-p
                                (genetics-kit-records-in kit chrom pos pos))))
    (or (and ref (cl-find ref hits :key #'genetics-snp-ref :test #'equal))
        (car hits))))

(defconst genetics--max-deletion 1000
  "Longest upstream deletion checked when inferring a reference call.")

(defun genetics-kit-site-spanned-p (kit chrom pos)
  "Return the non-reference record of KIT whose REF covers CHROM:POS.
Only records starting before POS (a deletion that spans it) count."
  (cl-find-if
   (lambda (snp)
     (let ((ref (genetics-snp-ref snp)))
       (and ref (< (genetics-snp-pos snp) pos)
            (>= (+ (genetics-snp-pos snp) (length ref) -1) pos)
            (not (equal (genetics-snp-genotype snp) (concat ref ref))))))
   (genetics-kit-records-in kit chrom (max 1 (- pos genetics--max-deletion))
                            (1- pos))))

;;;; Kit registry

(defun genetics--unique-name (name &optional file)
  "Return NAME, suffixed with <N> if a loaded kit already has it.
A kit loaded from FILE is being replaced and does not count."
  (let ((n 1) (candidate name)
        (others (cl-remove file genetics-loaded-kits
                           :key #'genetics-kit-file :test #'equal)))
    (while (cl-find candidate others :key #'genetics-kit-name :test #'equal)
      (setq n (1+ n) candidate (format "%s<%d>" name n)))
    candidate))

(defun genetics-register-kit (kit)
  "Add KIT to `genetics-loaded-kits' (replacing one for the same file)."
  (setq genetics-loaded-kits
        (cons kit (cl-remove (genetics-kit-file kit) genetics-loaded-kits
                             :key #'genetics-kit-file :test #'equal)))
  kit)

(defun genetics-find-kit (kit-or-name)
  "Return the loaded kit designated by KIT-OR-NAME (a kit or a name)."
  (cond ((genetics-kit-p kit-or-name) kit-or-name)
        ((cl-find kit-or-name genetics-loaded-kits
                  :key #'genetics-kit-name :test #'equal))
        (t (genetics--error 'genetics-no-kit "No loaded kit named %S"
                            kit-or-name))))

(defun genetics--read-kit (&optional prompt)
  "Return the kit for the current buffer, or ask using PROMPT."
  (cond (genetics--buffer-kit)
        ((null genetics-loaded-kits)
         (genetics--error 'genetics-no-kit
                          "No kit loaded; run M-x genetics-open first"))
        ((null (cdr genetics-loaded-kits)) (car genetics-loaded-kits))
        (t (genetics-find-kit
            (completing-read (or prompt "Kit: ")
                             (mapcar #'genetics-kit-name genetics-loaded-kits)
                             nil t nil nil
                             (genetics-kit-name (car genetics-loaded-kits)))))))

(defun genetics--format-label (format)
  "Return a human readable label for FORMAT symbol."
  (pcase format
    ('23andme "23andMe") ('ancestry "AncestryDNA")
    ('myheritage "MyHeritage CSV") ('ftdna "FamilyTreeDNA CSV")
    ('vcf "VCF") ('gvcf "gVCF") ('fastq-derived "VCF called from FASTQ")
    (_ (format "%s" format))))

(provide 'genetics-core)
;;; genetics-core.el ends here
