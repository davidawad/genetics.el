;;; genetics-detect.el --- Format, build and reference-call detection for genetics.el -*- lexical-binding: t; -*-

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

;; Recognizes 23andMe, AncestryDNA, MyHeritage/FTDNA CSV, VCF, FASTQ and
;; aligned-read files, the genome build, gVCFs and the VCF reference-call
;; model, and writes the strand note shown for each kit.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)

(defconst genetics--detect-bytes 65536
  "Number of leading bytes inspected when detecting a format.")


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

(defconst genetics--fastq-name-regexp
  "\\.\\(fastq\\|fq\\)\\(\\.gz\\|\\.bz2\\|\\.xz\\|\\.zst\\)?\\'"
  "File names that are FASTQ reads.")

(defconst genetics--alignment-name-regexp
  "\\.\\(bam\\|cram\\|sam\\)\\'"
  "File names that are aligned reads.")

(defun genetics-fastq-file-p (file)
  "Return non-nil if FILE is a FASTQ file, judged by name or content.
A FASTQ record is four lines: @name, bases, +, qualities."
  (or (let ((case-fold-search t))
        (string-match-p genetics--fastq-name-regexp file))
      (and (not (genetics--gz-file-p file))
           (file-readable-p file)
           (not (file-directory-p file))
           (let ((lines (genetics--head-lines
                         (with-temp-buffer
                           (set-buffer-multibyte nil)
                           (insert-file-contents-literally file nil 0 4096)
                           (buffer-string)))))
             (and (string-prefix-p "@" (or (nth 0 lines) ""))
                  (string-match-p "\\`[ACGTNacgtn.]+\\'" (or (nth 1 lines) ""))
                  (string-prefix-p "+" (or (nth 2 lines) "")))))))

(defun genetics--read-head (file)
  "Return the first `genetics--detect-bytes' bytes of FILE as a string."
  (with-temp-buffer
    (set-buffer-multibyte nil)
    (insert-file-contents-literally file nil 0 genetics--detect-bytes)
    (buffer-string)))

(defun genetics-detect-format (file)
  "Return the format symbol of FILE.
One of `23andme', `ancestry', `myheritage', `ftdna', `vcf', or `fastq'
and `alignment' for raw or aligned reads, which hold no genotypes.
Other compressed (.gz) files are assumed to be VCF and verified when
parsed."
  (unless (file-readable-p file)
    (genetics--error 'genetics-file-error "Cannot read file: %s" file))
  (cond ((genetics-fastq-file-p file) 'fastq)
        ((let ((case-fold-search t))
           (string-match-p genetics--alignment-name-regexp file))
         'alignment)
        ((genetics--gz-file-p file) 'vcf)
        (t (genetics--detect-head (genetics--read-head file)))))

(defconst genetics-fastq-explanation
  "FASTQ files hold raw sequencer reads, not genotypes. Reads must be aligned to a reference genome and variant-called to produce a VCF before they can be read here. If your provider (e.g. Nucleus) also gave you a VCF, open that: it is the product of these reads. Otherwise `genetics-fastq-plan' shows how genome-cli would do it (`genome pipeline plan') and `genetics-fastq-run' runs it."
  "Why a FASTQ file cannot be opened directly.")

(defun genetics--signal-unsupported (file format)
  "Signal an error: FILE of FORMAT has no genotypes, and say why."
  (pcase format
    ('fastq (genetics--error 'genetics-fastq-file "%s: %s"
                             (file-name-nondirectory file)
                             genetics-fastq-explanation))
    (_ (genetics--error 'genetics-unsupported-file
                        "%s holds aligned reads (BAM/CRAM/SAM), not genotypes; variant-call it to a VCF first"
                        (file-name-nondirectory file)))))

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

(defun genetics--gvcf-header-p (text)
  "Return non-nil if VCF header TEXT declares a gVCF (reference blocks)."
  (or (string-search "##GVCFBlock" text)
      (string-search "<NON_REF>" text)
      (string-search "ID=*,Description=\"Represents any possible" text)))

(defun genetics--vcf-ref-calls (gvcf count homref)
  "Decide the reference-call model of a VCF.
GVCF is non-nil for a gVCF header, COUNT the record count and HOMREF the
number of explicit homozygous-reference calls.  See
`genetics-vcf-ref-calls'.  Not cached: it follows the current options."
  (cond (gvcf 'explicit)
        ((not (eq genetics-vcf-ref-calls 'auto)) genetics-vcf-ref-calls)
        ((and (zerop homref) (>= count genetics-wgs-min-records))
         'absent-means-ref)
        (t 'unknown)))

(defun genetics--hom-ref-gt-p (line sample-index)
  "Return non-nil if VCF LINE has a homozygous-reference GT in SAMPLE-INDEX."
  (let* ((f (split-string (string-trim-right line "\r") "\t"))
         (fmt (nth 8 f)) (sample (nth sample-index f))
         (i (and fmt sample
                 (cl-position "GT" (split-string fmt ":") :test #'equal)))
         (gt (and i (nth i (split-string sample ":"))))
         (idx (and gt (split-string gt "[/|]" t))))
    (and idx (cl-every (lambda (a) (equal a "0")) idx))))

(defun genetics--apply-vcf-model (kit)
  "Set the assay and reference-call model of eager KIT from its stats."
  (if (eq (genetics-kit-format kit) 'vcf)
      (let* ((s (genetics-kit-stats kit))
             (rc (genetics--vcf-ref-calls (plist-get s :gvcf)
                                          (plist-get s :total)
                                          (plist-get s :hom-ref))))
        (setf (genetics-kit-ref-calls kit) rc
              (genetics-kit-assay kit)
              (genetics--vcf-assay rc (plist-get s :total))))
    (setf (genetics-kit-ref-calls kit) 'explicit
          (genetics-kit-assay kit) 'array))
  kit)

(defun genetics--vcf-assay (ref-calls count)
  "Guess the assay of a VCF from REF-CALLS and record COUNT."
  (if (or (eq ref-calls 'absent-means-ref)
          (>= count genetics-wgs-min-records))
      'wgs
    'unknown))

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

(provide 'genetics-detect)
;;; genetics-detect.el ends here
