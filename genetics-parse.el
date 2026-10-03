;;; genetics-parse.el --- Format detection and parsers for genetics.el -*- lexical-binding: t; -*-

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

;; Detects and parses 23andMe, AncestryDNA, MyHeritage/FTDNA CSV and VCF
;; (plain or .vcf.gz) files into a `genetics-kit'.  Large VCFs are
;; offset-indexed instead of loaded.  Compressed files are decompressed by
;; the local gzip executable; nothing is sent anywhere.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)
(require 'genetics-stats)

(defconst genetics--detect-bytes 65536
  "Number of leading bytes inspected when detecting a format.")

(defconst genetics--cache-version 1
  "Version stamp of the on-disk cache layout.")

;;;; Detection

(defun genetics--head-lines (head)
  "Split string HEAD into lines without CR characters."
  (split-string (replace-regexp-in-string "\r" "" head) "\n"))

(defun genetics--detect-head (head)
  "Return the format symbol for the leading text HEAD, or signal an error."
  (let* ((case-fold-search t)
         (lines (genetics--head-lines head))
         (comments (cl-remove-if-not (lambda (l) (string-prefix-p "#" l))
                                     lines))
         (data (cl-find-if (lambda (l) (and (not (string-empty-p l))
                                            (not (string-prefix-p "#" l))))
                           lines))
         (ctext (mapconcat #'identity comments "\n")))
    (cond
     ((string-prefix-p "##fileformat=VCF" head) 'vcf)
     ((null data)
      (genetics--error 'genetics-unknown-format
                       "No data lines found; not a supported genotype file"))
     ((and (not (string-search "\t" data)) (string-search "," data))
      (if (string-match-p "myheritage" ctext) 'myheritage 'ftdna))
     ((or (string-match-p "ancestrydna" ctext)
          (string-match-p "\\`rsid\t+chromosome\t+position\t+allele1" data)
          (= (length (split-string data "\t")) 5))
      'ancestry)
     ((or (string-match-p "23andme" ctext)
          (string-match-p "rsid\tchromosome\tposition\tgenotype" ctext)
          (= (length (split-string data "\t")) 4))
      '23andme)
     (t (genetics--error
         'genetics-unknown-format
         "Unrecognized file; expected 23andMe, AncestryDNA, CSV or VCF")))))

(defun genetics--gz-file-p (file)
  "Return non-nil if FILE has a .gz extension."
  (string-suffix-p ".gz" file t))

(defun genetics--read-head (file)
  "Return the first `genetics--detect-bytes' bytes of FILE as a string."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (insert-file-contents-literally file nil 0 genetics--detect-bytes)
    (buffer-string)))

(defun genetics-detect-format (file)
  "Return the format symbol of FILE.
One of `23andme', `ancestry', `myheritage', `ftdna' or `vcf'.  Compressed
\(.gz) files are assumed to be VCF and verified when parsed."
  (unless (file-readable-p file)
    (genetics--error 'genetics-file-error "Cannot read file: %s" file))
  (if (genetics--gz-file-p file)
      'vcf
    (genetics--detect-head (genetics--read-head file))))

(defun genetics--detect-build (text)
  "Return \"37\", \"38\" or nil from header TEXT."
  (let ((case-fold-search t))
    (cond ((string-match-p "grch38\\|hg38\\|build[ _]*38\\|\\bb38\\b" text)
           "38")
          ((string-match-p "grch37\\|hg19\\|build[ _]*37\\|hs37\\|\\bb37\\b"
                           text)
           "37")
          ((string-search "length=249250621" text) "37")
          ((string-search "length=248956422" text) "38"))))

(defun genetics--detect-chip (text count hint)
  "Return a chip version string (v3, v4, v5 or unknown).
TEXT is the header, COUNT the number of SNPs, HINT an explicit override."
  (cond (hint (format "%s" hint))
        ((let ((case-fold-search t))
           (and (string-match "\\bv\\([345]\\)\\b" text)
                (concat "v" (match-string 1 text)))))
        ((<= 610000 count 700000) "v5")
        ((and (<= 540000 count) (< count 610000)) "v4")
        ((<= 900000 count 1100000) "v3")
        (t "unknown")))

;;;; Strand notes

(defun genetics--strand-note (format build assumed)
  "Return the strand/build note for FORMAT and BUILD.
ASSUMED non-nil means BUILD was not declared in the file."
  (concat
   (pcase format
     ('23andme "23andMe reports genotypes on the + (forward) strand of GRCh37; no liftover is performed.")
     ('ancestry "AncestryDNA genotypes are reported on the + (forward) strand of build 37.1 (GRCh37); no liftover is performed.")
     ('myheritage "MyHeritage genotypes are reported on build 37 (GRCh37), assumed + strand; no liftover is performed.")
     ('ftdna "FamilyTreeDNA genotypes are reported on build 37 (GRCh37), assumed + strand; no liftover is performed.")
     ('vcf (format "VCF genotypes are expressed relative to the file's REF/ALT alleles on build %s; no liftover is performed."
                   (or build "unknown")))
     (_ ""))
   (if assumed " The build was not declared in the file; the vendor default was assumed." "")))

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
                 (unless (gethash (genetics-snp-rsid snp) table)
                   (puthash (genetics-snp-rsid snp) snp table)))))
        (forward-line 1)))
    (unless info
      (setq info (genetics--vcf-header-info (reverse header) sample)))
    (list (nreverse header) (nth 1 info) table (car info))))

(defun genetics--index-vcf (data-file sample)
  "Offset-index VCF DATA-FILE for SAMPLE without loading its records.
Return a plist with :index :ranges :counts :header :sample-index."
  (let ((index (make-hash-table :test 'equal))
        (ranges (make-hash-table :test 'equal))
        (counts (make-hash-table :test 'equal))
        (header nil))
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
             (puthash chrom (1+ (gethash chrom counts 0)) counts)
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
              :build (nth 0 info) :sample (nth 1 info)
              :sample-index (nth 2 info)
              :counts (mapcar (lambda (c) (cons c (gethash c counts)))
                              (genetics--sort-chroms
                               (hash-table-keys counts))))))))

;;;; gzip

(defun genetics--gzip ()
  "Return the gzip executable or signal `genetics-gzip-error'."
  (or (executable-find "gzip")
      (genetics--error 'genetics-gzip-error
                       "gzip executable not found; install gzip to read .gz files")))

(defun genetics--gunzip-into-buffer (file)
  "Insert the decompressed contents of FILE at point."
  (let* ((coding-system-for-read 'no-conversion)
         (status (call-process (genetics--gzip) nil t nil "-dc" "--" file)))
    (unless (eq status 0)
      (genetics--error 'genetics-gzip-error
                       "gzip failed on %s (exit status %s)" file status))))

(defun genetics--cache-key (file)
  "Return a hash string identifying FILE by truename, size and mtime."
  (let* ((truename (file-truename file))
         (attrs (file-attributes truename)))
    (secure-hash 'sha1 (format "%s|%d|%s" truename
                               (file-attribute-size attrs)
                               (format-time-string
                                "%s.%N" (file-attribute-modification-time
                                         attrs))))))

(defun genetics--decompressed-copy (file)
  "Return the path of a decompressed copy of gz FILE in the cache directory."
  (let ((out (expand-file-name (format "vcf-%s.vcf" (genetics--cache-key file))
                               genetics-cache-directory)))
    (unless (and (file-exists-p out)
                 (> (file-attribute-size (file-attributes out)) 0))
      (make-directory genetics-cache-directory t)
      (let ((tmp (concat out ".part")))
        (let ((status (call-process (genetics--gzip) nil (list :file tmp) nil
                                    "-dc" "--" file)))
          (unless (eq status 0)
            (ignore-errors (delete-file tmp))
            (genetics--error 'genetics-gzip-error
                             "gzip failed on %s (exit status %s)" file status)))
        (rename-file tmp out t)))
    out))

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
    (with-temp-file (genetics--cache-file key)
      (let ((print-length nil) (print-level nil) (coding-system-for-write 'utf-8))
        (prin1 (list 'genetics-cache genetics--cache-version
                     :key (genetics--cache-key key)
                     :format (genetics-kit-format kit)
                     :build (genetics-kit-build kit)
                     :chip (genetics-kit-chip kit)
                     :sample (genetics-kit-sample kit)
                     :note (genetics-kit-strand-note kit)
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
              :strand-note (plist-get p :note) :table table))))))))

;;;; Entry point

(defun genetics--finish-kit (kit)
  "Index and compute statistics for eager KIT, then return it."
  (setf (genetics-kit-chroms kit)
        (genetics--sort-chroms-into (genetics-kit-table kit)))
  (setf (genetics-kit-stats kit) (genetics-compute-stats kit))
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
             (build (or build (and (not (eq format 'vcf)) "37"))))
        (genetics--finish-kit
         (genetics-kit--create
          :file file :format format :build build
          :chip (if (eq format '23andme)
                    (genetics--detect-chip header (hash-table-count table)
                                           chip-hint)
                  "unknown")
          :sample sname
          :strand-note (genetics--strand-note format build assumed)
          :table table))))))

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
               :line-parser (lambda (line) (genetics--parse-vcf-line line sidx)))))
    (unless (car (plist-get info :ranges))
      (genetics--error 'genetics-parse-error "No genotype records in %s" file))
    (setf (genetics-kit-stats kit)
          (list :total (hash-table-count (plist-get info :index))
                :chrom-counts counts :lazy t))
    kit))

;;;###autoload
(defun genetics-parse-file (file &rest args)
  "Parse genotype FILE and return a `genetics-kit' (not registered).
ARGS is a plist: :format forces a format symbol, :chip gives a chip version
hint, :sample selects a VCF sample column, :lazy non-nil forces offset-index
mode for VCF, :cache nil disables the on-disk cache for this call."
  (let* ((file (expand-file-name file))
         (format (or (plist-get args :format) (genetics-detect-format file)))
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
    kit))

(provide 'genetics-parse)
;;; genetics-parse.el ends here
