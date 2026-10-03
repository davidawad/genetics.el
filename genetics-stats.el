;;; genetics-stats.el --- Statistics, sex inference and caveats for genetics.el -*- lexical-binding: t; -*-

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

;; Summary statistics over a kit (no-calls, heterozygosity, per-chromosome
;; counts), chromosomal sex inference, and the strand/build caveat text.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)

(defun genetics-compute-stats (kit)
  "Compute the statistics plist for eager KIT.
Keys: :total :nocalls :het :hom :hemi :other :chrom-counts and the sex
marker counts :x-het :x-hom :x-hemi :y-total :y-called."
  (let ((nocalls 0) (het 0) (hom 0) (hemi 0) (other 0)
        (xhet 0) (xhom 0) (xhemi 0) (ytotal 0) (ycalled 0)
        (counts (make-hash-table :test 'eq)))
    (maphash
     (lambda (_k snp)
       (let ((chrom (genetics-snp-chrom snp))
             (z (genetics-zygosity (genetics-snp-genotype snp))))
         (puthash chrom (1+ (gethash chrom counts 0)) counts)
         (pcase z
           ('no-call (cl-incf nocalls))
           ('heterozygous (cl-incf het))
           ('homozygous (cl-incf hom))
           ('hemizygous (cl-incf hemi))
           (_ (cl-incf other)))
         (cond ((equal chrom "X")
                (pcase z
                  ('heterozygous (cl-incf xhet))
                  ('homozygous (cl-incf xhom))
                  ('hemizygous (cl-incf xhemi))))
               ((equal chrom "Y")
                (cl-incf ytotal)
                (unless (eq z 'no-call) (cl-incf ycalled))))))
     (genetics-kit-table kit))
    (list :total (hash-table-count (genetics-kit-table kit))
          :nocalls nocalls :het het :hom hom :hemi hemi :other other
          :chrom-counts (mapcar (lambda (c) (cons c (gethash c counts)))
                                (genetics--sort-chroms
                                 (hash-table-keys counts)))
          :x-het xhet :x-hom xhom :x-hemi xhemi
          :y-total ytotal :y-called ycalled)))

(defun genetics-nocall-rate (kit)
  "Return the no-call fraction of KIT, or nil when unknown."
  (let ((s (genetics-kit-stats kit)))
    (when (and (plist-get s :nocalls) (> (plist-get s :total) 0))
      (/ (float (plist-get s :nocalls)) (plist-get s :total)))))

(defun genetics-sex-evidence (kit)
  "Return a plist of X heterozygosity and Y call rates for KIT."
  (let* ((s (genetics-kit-stats kit))
         (xhet (or (plist-get s :x-het) 0))
         (xhom (or (plist-get s :x-hom) 0))
         (xhemi (or (plist-get s :x-hemi) 0))
         (diploid (+ xhet xhom))
         (ytotal (or (plist-get s :y-total) 0)))
    (list :x-het-rate (cond ((> diploid 0) (/ (float xhet) diploid))
                            ((> xhemi 0) 0.0))
          :x-calls (+ diploid xhemi)
          :y-called-rate (when (> ytotal 0)
                           (/ (float (plist-get s :y-called)) ytotal)))))

(defun genetics-infer-sex (kit)
  "Infer chromosomal gender of KIT: `male', `female' or `uncertain'.
Female: X heterozygosity above 1% and mostly no-called Y.  Male: Y calls
present and X heterozygosity below 1%.  Offset-indexed kits are uncertain."
  (if (genetics-kit-lazy kit)
      'uncertain
    (let* ((e (genetics-sex-evidence kit))
           (xr (plist-get e :x-het-rate))
           (yr (plist-get e :y-called-rate)))
      (cond ((null xr) 'uncertain)
            ((and (> xr 0.01) (or (null yr) (< yr 0.2))) 'female)
            ((and (< xr 0.01) yr (>= yr 0.5)) 'male)
            (t 'uncertain)))))

(defun genetics-sex-description (kit)
  "Return a sentence describing the gender inference for KIT."
  (let* ((e (genetics-sex-evidence kit))
         (xr (plist-get e :x-het-rate))
         (yr (plist-get e :y-called-rate)))
    (format "%s (X heterozygosity %s, Y called %s)"
            (genetics-infer-sex kit)
            (if xr (format "%.2f%%" (* 100 xr)) "n/a")
            (if yr (format "%.0f%%" (* 100 yr)) "n/a"))))

(defun genetics-kit-caveats (kit)
  "Return a list of caveat strings (strand, build, chip) for KIT."
  (delq nil
        (list
         (genetics-kit-strand-note kit)
         (format "Positions are on GRCh%s; they cannot be compared to data on another build without liftover, which this package does not do."
                 (or (genetics-kit-build kit) "(unknown build)"))
         (when (equal (genetics-kit-chip kit) "unknown")
           (when (eq (genetics-kit-format kit) '23andme)
             "The 23andMe chip version could not be determined."))
         "Palindromic SNPs (A/T and C/G) cannot have their strand verified from the genotype alone.")))

(defun genetics-builds-differ-p (a b)
  "Return non-nil if kits A and B have different known builds."
  (let ((ba (genetics-kit-build a)) (bb (genetics-kit-build b)))
    (and ba bb (not (equal ba bb)))))

(provide 'genetics-stats)
;;; genetics-stats.el ends here
