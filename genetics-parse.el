;;; genetics-parse.el --- Format detection and parsers for genetics.el -*- lexical-binding: t; -*-

;; Copyright (C) 2026 David Awad

;; Author: David Awad <davidawad@protonmail.com>
;; Maintainer: David Awad <davidawad@protonmail.com>
;; Version: 0.1.0
;; Package-Requires: ((emacs "29.1"))
;; Keywords: data, tools
;; URL: https://github.com/davidawad/genetics.el

;; This file is not part of GNU Emacs.

;; Permission is hereby granted, free of charge, to any person obtaining a
;; copy of this software, to deal in it without restriction (MIT licence).

;;; Commentary:

;; Detects and parses 23andMe, AncestryDNA, MyHeritage/FTDNA CSV and VCF
;; (plain or .vcf.gz) files into a `genetics-kit'.  Large VCFs are
;; offset-indexed instead of loaded.  Compressed files are decompressed by
;; the local gzip executable (or Emacs' zlib); nothing is sent anywhere.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)
(require 'genetics-stats)
(require 'genetics-detect)
(require 'genetics-gzip)

(defconst genetics--block-size 256
  "Records between two entries of the position index of a lazy VCF.")

(defconst genetics--cache-version 3
  "Version stamp of the on-disk cache layout.")

;;;; Field scanning

(defsubst genetics--field (eol sep)
  "Return the field at point up to SEP or EOL and move past the separator.
Surrounding double quotes are removed."
  (let ((s (point)))
    (skip-chars-forward (if (eq sep ?\t) "^\t" "^,") eol)
    (let ((str (buffer-substring-no-properties s (point))))
      (when (< (point) eol) (forward-char 1))
      (if (and (> (length str) 1) (eq (aref str 0) ?\"))
          (substring str 1 -1)
        str))))

(defun genetics--norm-genotype (s)
  "Normalize vendor genotype string S (no-call becomes `--')."
  (if (member s '("" "--" "0" "00" "-" "NC"))
      "--"
    (upcase s)))

(defun genetics--sort-chroms-into (table)
  "Return an alist (CHROM . sorted vector of SNPs) built from hash TABLE."
  (let ((groups (make-hash-table :test 'eq)))
    (maphash (lambda (_k snp) (push snp (gethash (genetics-snp-chrom snp)
                                                 groups)))
             table)
    (let (result)
      (dolist (c (genetics--sort-chroms (hash-table-keys groups)))
        (push (cons c (genetics--sort-vector (gethash c groups))) result))
      (nreverse result))))

(defun genetics--sort-vector (snps)
  "Return a vector of SNPS (list) sorted by position then rsid."
  (let ((vec (vconcat snps)))
    (sort vec (lambda (a b)
                (let ((pa (genetics-snp-pos a)) (pb (genetics-snp-pos b)))
                  (if (= pa pb)
                      (string< (genetics-snp-rsid a) (genetics-snp-rsid b))
                    (< pa pb)))))))

;;;; Tab / CSV parsing (23andMe, AncestryDNA, MyHeritage, FTDNA)

(defun genetics--parse-array-buffer (format)
  "Parse the current buffer (FORMAT: 23andme ancestry myheritage ftdna).
Return (HEADER-TEXT . HASH-TABLE of rsid to SNP)."
  (goto-char (point-min))
  (let* ((csv (memq format '(myheritage ftdna)))
         (sep (if csv ?, ?\t))
         (ancestry (eq format 'ancestry))
         (table (make-hash-table :test 'equal :size 4096))
         (header nil) (seen-data nil) (lineno 0))
    (while (not (eobp))
      (setq lineno (1+ lineno))
      (let ((eol (line-end-position)) (c (following-char)))
        (when (and (> eol (point)) (eq (char-before eol) ?\r))
          (setq eol (1- eol)))
        (cond
         ((= (point) eol))
         ((= c ?#)
          (push (buffer-substring-no-properties (point) eol) header))
         (t
          (let* ((rsid (genetics--field eol sep))
                 (chrom (genetics--field eol sep))
                 (pos (genetics--field eol sep))
                 (g1 (genetics--field eol sep))
                 (g2 (when ancestry (genetics--field eol sep))))
            (cond
             ((and (not seen-data) (> (length pos) 0)
                   (not (<= ?0 (aref pos 0) ?9)))
              (setq seen-data t))       ; column header line
             ((or (zerop (length pos)) (not (<= ?0 (aref pos 0) ?9))
                  (zerop (length chrom)) (zerop (length rsid)))
              (genetics--error 'genetics-parse-error
                               "Malformed line %d (need rsid, chromosome, position, genotype)"
                               lineno))
             (t
              (setq seen-data t)
              (puthash rsid
                       (genetics-snp-create
                        rsid (genetics-normalize-chrom chrom ancestry)
                        (string-to-number pos)
                        (if ancestry
                            (if (or (member g1 '("" "0")) (member g2 '("" "0")))
                                "--"
                              (upcase (concat g1 g2)))
                          (genetics--norm-genotype g1)))
                       table)))))))
      (forward-line 1))
    (cons (mapconcat #'identity (nreverse header) "\n") table)))

;;;; VCF

(defun genetics--gt-genotype (gt ref alts)
  "Convert VCF genotype GT to a genotype string using REF and ALTS (list)."
  (let ((idx (split-string gt "[/|]" t)))
    (if (or (null idx) (member "." idx))
        "--"
      (let ((letters
             (mapcar
              (lambda (i)
                (unless (string-match-p "\\`[0-9]+\\'" i)
                  (genetics--error 'genetics-parse-error
                                   "Bad GT value %S" gt))
                (let ((n (string-to-number i)))
                  (if (= n 0)
                      ref
                    (or (nth (1- n) alts)
                        (genetics--error 'genetics-parse-error
                                         "GT %S refers to a missing ALT allele"
                                         gt)))))
              idx)))
        (upcase (if (cl-every (lambda (s) (= (length s) 1)) letters)
                    (apply #'concat letters)
                  (mapconcat #'identity letters "/")))))))

(defun genetics--parse-vcf-line (line sample-index)
  "Parse VCF data LINE using the genotype column SAMPLE-INDEX.
Return a `genetics-snp'.  Non-rsid ids are stored as chrom:pos."
  (when (string-suffix-p "\r" line)
    (setq line (substring line 0 -1)))
  (let ((f (split-string line "\t")))
    (when (< (length f) 8)
      (genetics--error 'genetics-parse-error "Malformed VCF line: %.60s" line))
    (let* ((chrom (genetics-normalize-chrom (nth 0 f)))
           (pos (string-to-number (nth 1 f)))
           (id (nth 2 f))
           (ref (upcase (nth 3 f)))
           (alts (unless (member (nth 4 f) '("." "")) (split-string (nth 4 f) ",")))
           (fmt (nth 8 f))
           (sample (nth sample-index f))
           (key (if (member id '("." "")) (format "%s:%d" chrom pos)
                  (car (split-string id ";"))))
           (gt (when (and fmt sample)
                 (let ((i (cl-position "GT" (split-string fmt ":") :test #'equal)))
                   (when i (nth i (split-string sample ":")))))))
      (genetics-snp-create key chrom pos
                           (if gt (genetics--gt-genotype gt ref alts) "--")
                           ref (when alts (mapconcat #'identity alts ","))))))

(defun genetics--vcf-header-info (lines &optional sample)
  "Return (BUILD SAMPLE-NAME SAMPLE-INDEX) from VCF header LINES.
SAMPLE selects a sample column by name; default is the first."
  (let* ((text (mapconcat #'identity lines "\n"))
         (chrom-line (cl-find-if (lambda (l) (string-prefix-p "#CHROM" l))
                                 lines))
         (cols (and chrom-line (split-string (string-trim-right chrom-line)
                                             "\t")))
         (names (nthcdr 9 cols))
         (name (if sample
                   (or (car (member sample names))
                       (genetics--error 'genetics-parse-error
                                        "VCF has no sample %S (have: %s)"
                                        sample (string-join names ", ")))
                 (car names))))
    (unless chrom-line
      (genetics--error 'genetics-parse-error "VCF has no #CHROM header line"))
    (unless name
      (genetics--error 'genetics-parse-error
                       "VCF has no sample (genotype) column"))
    (list (genetics--detect-build text) name
          (+ 9 (cl-position name names :test #'equal)))))

(defun genetics--vcf-table-add (table snp)
  "Add VCF record SNP to TABLE under its id without dropping records.
A second record with the same id (e.g. an indel and a SNP at one
position) is kept under id~REF>ALT; the single-base record keeps the
plain id so lookups by rsid or chrom:pos find the SNP."
  (let* ((key (genetics-snp-rsid snp))
         (old (gethash key table)))
    (if (null old)
        (puthash key snp table)
      (when (and (genetics--snv-record-p snp)
                 (not (genetics--snv-record-p old)))
        (puthash key snp table)
        (setq snp old))
      (let* ((base (format "%s~%s>%s" key (genetics-snp-ref snp)
                           (genetics-snp-alt snp)))
             (alt base) (n 1))
        (while (gethash alt table)
          (setq n (1+ n) alt (format "%s#%d" base n)))
        (setf (genetics-snp-rsid snp) alt)
        (puthash alt snp table)))))

(defun genetics--parse-vcf-buffer (sample)
  "Parse the VCF in the current buffer for SAMPLE.
Return (HEADER-LINES SAMPLE-NAME TABLE BUILD)."
  (goto-char (point-min))
  (let ((header nil) (info nil) (table (make-hash-table :test 'equal :size 4096)))
    (while (not (eobp))
      (let ((eol (line-end-position)) (c (following-char)))
        (cond ((= (point) eol))
              ((= c ?#)
               (push (buffer-substring-no-properties (point) eol) header))
              (t
               (unless info
                 (setq info (genetics--vcf-header-info (reverse header) sample)))
               (let ((snp (genetics--parse-vcf-line
                           (buffer-substring-no-properties (point) eol)
                           (nth 2 info))))
                 (genetics--vcf-table-add table snp))))
        (forward-line 1)))
    (unless info
      (setq info (genetics--vcf-header-info (reverse header) sample)))
    (list (nreverse header) (nth 1 info) table (car info))))

(defun genetics--index-vcf (data-file sample)
  "Offset-index VCF DATA-FILE for SAMPLE without loading its records.
Return a plist with :index :ranges :counts :header :sample-index
:hom-ref :rsids and :blocks, an alist (CHROM . vector of (POS . OFFSET))
holding every `genetics--block-size'th record, used to seek by position."
  (let ((index (make-hash-table :test 'equal))
        (ranges (make-hash-table :test 'equal))
        (counts (make-hash-table :test 'equal))
        (blocks (make-hash-table :test 'equal))
        (header nil) (homref 0) (rsids 0) (sidx nil))
    (genetics--stream-lines
     data-file 0 nil
     (lambda (line off)
       (cond
        ((string-empty-p line))
        ((eq (aref line 0) ?#) (push line header))
        (t
         (let* ((t1 (string-search "\t" line))
                (t2 (and t1 (string-search "\t" line (1+ t1))))
                (t3 (and t2 (string-search "\t" line (1+ t2)))))
           (unless t3
             (genetics--error 'genetics-parse-error
                              "Malformed VCF line at byte %d" off))
           (let* ((chrom (genetics-normalize-chrom (substring line 0 t1)))
                  (pos (substring line (1+ t1) t2))
                  (id (substring line (1+ t2) t3))
                  (end (+ off (length line) 1))
                  (range (gethash chrom ranges)))
             (if range
                 (setcdr range end)
               (puthash chrom (cons off end) ranges))
             (when (zerop (% (gethash chrom counts 0)
                             genetics--block-size))
               (push (cons (string-to-number pos) off)
                     (gethash chrom blocks)))
             (puthash chrom (1+ (gethash chrom counts 0)) counts)
             (unless sidx
               (setq sidx (nth 2 (genetics--vcf-header-info (reverse header)
                                                            sample))))
             (when (and (or (string-search "0/0" line t3)
                            (string-search "0|0" line t3))
                        (genetics--hom-ref-gt-p line sidx))
               (cl-incf homref))
             (when (string-prefix-p "rs" id) (cl-incf rsids))
             (if (member id '("." ""))
                 (puthash (format "%s:%s" chrom pos) off index)
               (dolist (one (split-string id ";" t))
                 (unless (gethash one index)
                   (puthash one off index))))))))))
    (let ((hl (nreverse header)) (sorted nil))
      (dolist (c (genetics--sort-chroms (hash-table-keys ranges)))
        (push (cons c (gethash c ranges)) sorted))
      (let ((info (genetics--vcf-header-info hl sample)))
        (list :index index :ranges (nreverse sorted) :header hl
              :hom-ref homref :rsids rsids
              :blocks (mapcar (lambda (c)
                                (cons c (vconcat (nreverse (gethash c blocks)))))
                              (hash-table-keys blocks))
              :build (nth 0 info) :sample (nth 1 info)
              :sample-index (nth 2 info)
              :counts (mapcar (lambda (c) (cons c (gethash c counts)))
                              (genetics--sort-chroms
                               (hash-table-keys counts))))))))

;;;; Cache

(defun genetics--cache-file (file)
  "Return the cache file path for FILE."
  (expand-file-name (concat (genetics--cache-key file) ".eld")
                    genetics-cache-directory))

(defun genetics--cache-write (kit key)
  "Write eager KIT to the cache under KEY (a source file)."
  (let ((records nil))
    (maphash (lambda (_k s)
               (push (vector (genetics-snp-rsid s) (genetics-snp-chrom s)
                             (genetics-snp-pos s) (genetics-snp-genotype s)
                             (genetics-snp-ref s) (genetics-snp-alt s))
                     records))
             (genetics-kit-table kit))
    (make-directory genetics-cache-directory t)
    (genetics--with-output-file (genetics--cache-file key)
      (let ((print-length nil) (print-level nil))
        (prin1 (list 'genetics-cache genetics--cache-version
                     :key (genetics--cache-key key)
                     :format (genetics-kit-format kit)
                     :build (genetics-kit-build kit)
                     :chip (genetics-kit-chip kit)
                     :sample (genetics-kit-sample kit)
                     :note (genetics-kit-strand-note kit)
                     :gvcf (plist-get (genetics-kit-stats kit) :gvcf)
                     :records (vconcat records))
               (current-buffer))))))

(defun genetics--cache-read (file)
  "Return a kit rebuilt from the cache for FILE, or nil if absent or stale."
  (let ((cache (genetics--cache-file file)))
    (when (file-readable-p cache)
      (let ((data (with-temp-buffer
                    (let ((coding-system-for-read 'utf-8))
                      (insert-file-contents cache))
                    (ignore-errors (read (current-buffer))))))
        (when (and (eq (car-safe data) 'genetics-cache)
                   (eql (cadr data) genetics--cache-version)
                   (equal (plist-get (cddr data) :key)
                          (genetics--cache-key file)))
          (let ((table (make-hash-table :test 'equal :size 4096))
                (p (cddr data)))
            (cl-loop for r across (plist-get p :records)
                     do (puthash (aref r 0)
                                 (genetics-snp-create
                                  (aref r 0) (genetics--intern-chrom (aref r 1))
                                  (aref r 2) (aref r 3) (aref r 4) (aref r 5))
                                 table))
            (genetics--finish-kit
             (genetics-kit--create
              :file file :format (plist-get p :format)
              :build (plist-get p :build) :chip (plist-get p :chip)
              :sample (plist-get p :sample)
              :strand-note (plist-get p :note) :table table)
             (plist-get p :gvcf))))))))

;;;; Entry point

(defun genetics--finish-kit (kit &optional gvcf)
  "Index and compute statistics for eager KIT, then return it.
GVCF non-nil records that the VCF header declared a gVCF."
  (setf (genetics-kit-chroms kit)
        (genetics--sort-chroms-into (genetics-kit-table kit)))
  (setf (genetics-kit-stats kit)
        (append (genetics-compute-stats kit) (list :gvcf gvcf)))
  kit)

(defun genetics--base-name (file)
  "Return a display name for FILE without directory or extensions."
  (let ((n (file-name-nondirectory file)))
    (replace-regexp-in-string "\\(\\.gz\\)?\\(\\.[A-Za-z0-9]+\\)?\\'" "" n)))

(defun genetics--parse-eager (file format gz chip-hint sample)
  "Parse FILE of FORMAT fully into memory.  GZ, CHIP-HINT, SAMPLE as usual."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (if gz
        (genetics--gunzip-into-buffer file)
      (insert-file-contents-literally file))
    (when (zerop (buffer-size))
      (genetics--error 'genetics-parse-error "File is empty: %s" file))
    (when (and gz (not (string-prefix-p "##fileformat=VCF"
                                        (buffer-substring-no-properties
                                         (point-min)
                                         (min (point-max) 40)))))
      (genetics--error 'genetics-unknown-format
                       "Compressed file is not a VCF: %s" file))
    (let (header table build sname)
      (if (eq format 'vcf)
          (pcase-let ((`(,hl ,sn ,tb ,bd) (genetics--parse-vcf-buffer sample)))
            (setq header (mapconcat #'identity hl "\n") table tb
                  build bd sname sn))
        (pcase-let ((`(,h . ,tb) (genetics--parse-array-buffer format)))
          (setq header h table tb build (genetics--detect-build h))))
      (when (zerop (hash-table-count table))
        (genetics--error 'genetics-parse-error "No genotype records in %s" file))
      (let* ((assumed (and (null build) (not (eq format 'vcf))))
             (build (or build (and (not (eq format 'vcf)) "37")))
             (kit (genetics--finish-kit
                   (genetics-kit--create
                    :file file :format format :build build
                    :chip (if (eq format '23andme)
                              (genetics--detect-chip
                               header (hash-table-count table) chip-hint)
                            "unknown")
                    :sample sname
                    :strand-note (genetics--strand-note format build assumed)
                    :table table)
                   (and (eq format 'vcf) (genetics--gvcf-header-p header)))))
        kit))))

(defun genetics--parse-lazy (file gz sample)
  "Offset-index VCF FILE (GZ non-nil if compressed) for SAMPLE."
  (let* ((data (if gz (genetics--decompressed-copy file) file))
         (info (genetics--index-vcf data sample))
         (sidx (plist-get info :sample-index))
         (counts (plist-get info :counts))
         (kit (genetics-kit--create
               :file file :format 'vcf :build (plist-get info :build)
               :chip "unknown" :sample (plist-get info :sample)
               :strand-note (genetics--strand-note
                             'vcf (plist-get info :build) nil)
               :lazy t :index (plist-get info :index)
               :ranges (plist-get info :ranges) :data-file data
               :blocks (plist-get info :blocks)
               :line-parser (lambda (line) (genetics--parse-vcf-line line sidx)))))
    (unless (car (plist-get info :ranges))
      (genetics--error 'genetics-parse-error "No genotype records in %s" file))
    (let ((total (hash-table-count (plist-get info :index)))
          (rc (genetics--vcf-ref-calls
               (genetics--gvcf-header-p
                (mapconcat #'identity (plist-get info :header) "\n"))
               (hash-table-count (plist-get info :index))
               (plist-get info :hom-ref))))
      (setf (genetics-kit-stats kit)
            (list :total total :chrom-counts counts :lazy t
                  :rsids (plist-get info :rsids))
            (genetics-kit-ref-calls kit) rc
            (genetics-kit-assay kit) (genetics--vcf-assay rc total)))
    kit))

;;;###autoload
(defun genetics-parse-file (file &rest args)
  "Parse genotype FILE and return a `genetics-kit' (not registered).
ARGS is a plist: :format forces a format symbol, :chip gives a chip version
hint, :sample selects a VCF sample column, :lazy non-nil forces offset-index
mode for VCF, :cache nil disables the on-disk cache for this call,
:ref-calls and :assay override the detected reference-call model and
assay (see `genetics-kit').  FASTQ and BAM/CRAM files signal
`genetics-unsupported-file' with an explanation."
  (let* ((file (expand-file-name file))
         (format (or (plist-get args :format) (genetics-detect-format file)))
         (_ (when (memq format '(fastq alignment))
              (genetics--signal-unsupported file format)))
         (gz (genetics--gz-file-p file))
         (size (file-attribute-size (file-attributes file)))
         (lazy (and (eq format 'vcf)
                    (or (plist-get args :lazy) (> size genetics-vcf-eager-limit))))
         (use-cache (and genetics-use-cache (not lazy)
                         (not (plist-member args :chip))
                         (not (plist-get args :sample))
                         (if (plist-member args :cache) (plist-get args :cache) t)))
         (kit (or (and use-cache (genetics--cache-read file))
                  (let ((k (if lazy
                               (genetics--parse-lazy file gz (plist-get args :sample))
                             (genetics--parse-eager file format gz
                                                    (plist-get args :chip)
                                                    (plist-get args :sample)))))
                    (when use-cache (genetics--cache-write k file))
                    k))))
    (setf (genetics-kit-file kit) file)
    (setf (genetics-kit-name kit) (genetics--base-name file))
    (unless (genetics-kit-lazy kit)
      (genetics--apply-vcf-model kit))
    (when (plist-get args :ref-calls)
      (setf (genetics-kit-ref-calls kit) (plist-get args :ref-calls))
      (when (and (eq (plist-get args :ref-calls) 'absent-means-ref)
                 (eq (genetics-kit-assay kit) 'unknown))
        (setf (genetics-kit-assay kit) 'wgs)))
    (when (plist-get args :assay)
      (setf (genetics-kit-assay kit) (plist-get args :assay)))
    kit))

(provide 'genetics-parse)
;;; genetics-parse.el ends here
