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
                                                &optional ref alt))
                            (:copier nil))
  "One genotype call.  GENOTYPE is upper case letters or `--' for no-call."
  rsid chrom pos genotype ref alt)

(cl-defstruct (genetics-kit (:constructor genetics-kit--create)
                            (:copier nil))
  "A loaded raw-data file."
  name file format build chip strand-note sample
  table chroms stats lazy index ranges data-file line-parser)

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
  (if (genetics-kit-lazy kit)
      (let ((off (gethash rsid (genetics-kit-index kit))))
        (when off
          (funcall (genetics-kit-line-parser kit)
                   (genetics--read-line-at (genetics-kit-data-file kit) off))))
    (gethash rsid (genetics-kit-table kit))))

(defun genetics-kit-ids (kit)
  "Return a hash table whose keys are the ids available in KIT."
  (if (genetics-kit-lazy kit)
      (genetics-kit-index kit)
    (genetics-kit-table kit)))

(defun genetics-kit-snp-count (kit)
  "Return the number of records in KIT."
  (hash-table-count (genetics-kit-ids kit)))

(defun genetics-kit-chromosomes (kit)
  "Return the chromosomes present in KIT in canonical order."
  (if (genetics-kit-lazy kit)
      (genetics--sort-chroms (mapcar #'car (genetics-kit-ranges kit)))
    (mapcar #'car (genetics-kit-chroms kit))))

(defun genetics-kit-map-snps (kit fn &optional chrom)
  "Call FN on each SNP of KIT in chromosome/position order.
Only CHROM is visited when given.  FN returning `stop' ends the walk."
  (let ((chroms (if chrom (list chrom) (genetics-kit-chromosomes kit))))
    (catch 'genetics--stop
      (dolist (c chroms)
        (if (genetics-kit-lazy kit)
            (genetics--lazy-map-chrom kit c fn)
          (let ((vec (cdr (assoc c (genetics-kit-chroms kit)))))
            (cl-loop for snp across vec
                     when (eq (funcall fn snp) 'stop)
                     do (throw 'genetics--stop nil))))))))

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
    ('vcf "VCF") (_ (format "%s" format))))

(provide 'genetics-core)
;;; genetics-core.el ends here
