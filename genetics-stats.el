;;; genetics-stats.el --- Statistics, sex inference and caveats for genetics.el -*- lexical-binding: t; -*-

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

;; Summary statistics over a kit (no-calls, heterozygosity, per-chromosome
;; counts), chromosomal sex inference, and the strand/build caveat text.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)

;;;; Pseudoautosomal regions

(defconst genetics-par-regions
  '(("37" (60001 . 2699520) (154931044 . 155260560))
    ("38" (10001 . 2781479) (155701383 . 156030895)))
  "Pseudoautosomal regions PAR1 and PAR2 of chromosome X, per build.
Inclusive (START . END) positions from the Genome Reference Consortium.
Males are diploid there, so these sites are excluded from X
heterozygosity when inferring sex.")

(defun genetics-par-p (pos build)
  "Return non-nil if X position POS lies in a pseudoautosomal region.
BUILD is \"37\" or \"38\"; for any other value the regions of both builds
are excluded, which only drops a few extra sites."
  (cl-some (lambda (entry)
             (and (or (null build) (not (assoc build genetics-par-regions))
                      (equal build (car entry)))
                  (cl-some (lambda (r) (<= (car r) pos (cdr r))) (cdr entry))))
           genetics-par-regions))

;;;; Statistics

(defun genetics-compute-stats (kit)
  "Compute the statistics plist for eager KIT.
Keys: :total :nocalls :het :hom :hemi :other :hom-ref :rsids
:chrom-counts and the sex marker counts :x-het :x-hom :x-hemi :x-par
:y-total :y-called.  X sites in the pseudoautosomal regions of the kit's
build are counted in :x-par only, not in the X heterozygosity counts."
  (let ((nocalls 0) (het 0) (hom 0) (hemi 0) (other 0) (homref 0) (rsids 0)
        (xhet 0) (xhom 0) (xhemi 0) (xpar 0) (ytotal 0) (ycalled 0)
        (build (genetics-kit-build kit))
        (counts (make-hash-table :test 'eq)))
    (maphash
     (lambda (_k snp)
       (let* ((chrom (genetics-snp-chrom snp))
              (gt (genetics-snp-genotype snp))
              (z (genetics-zygosity gt))
              (ref (genetics-snp-ref snp)))
         (puthash chrom (1+ (gethash chrom counts 0)) counts)
         (when (string-prefix-p "rs" (genetics-snp-rsid snp))
           (cl-incf rsids))
         (pcase z
           ('no-call (cl-incf nocalls))
           ('heterozygous (cl-incf het))
           ('homozygous (cl-incf hom)
                        (when (and ref (equal gt (concat ref ref)))
                          (cl-incf homref)))
           ('hemizygous (cl-incf hemi))
           (_ (cl-incf other)))
         (cond ((equal chrom "X")
                (if (genetics-par-p (genetics-snp-pos snp) build)
                    (cl-incf xpar)
                  (pcase z
                    ('heterozygous (cl-incf xhet))
                    ('homozygous (cl-incf xhom))
                    ('hemizygous (cl-incf xhemi)))))
               ((equal chrom "Y")
                (cl-incf ytotal)
                (unless (eq z 'no-call) (cl-incf ycalled))))))
     (genetics-kit-table kit))
    (list :total (hash-table-count (genetics-kit-table kit))
          :nocalls nocalls :het het :hom hom :hemi hemi :other other
          :hom-ref homref :rsids rsids
          :chrom-counts (mapcar (lambda (c) (cons c (gethash c counts)))
                                (genetics--sort-chroms
                                 (hash-table-keys counts)))
          :x-het xhet :x-hom xhom :x-hemi xhemi :x-par xpar
          :y-total ytotal :y-called ycalled)))

(defun genetics-nocall-rate (kit)
  "Return the no-call fraction of KIT, or nil when unknown."
  (let ((s (genetics-kit-stats kit)))
    (when (and (plist-get s :nocalls) (> (plist-get s :total) 0))
      (/ (float (plist-get s :nocalls)) (plist-get s :total)))))

;;;; Sex inference

(defun genetics-sex-evidence (kit)
  "Return a plist of X heterozygosity and Y call rates for KIT.
Keys: :x-het-rate (outside the pseudoautosomal regions) :x-calls
:y-called-rate :y-called :method.  Kits from genome-cli carry the
values computed there."
  (let ((s (genetics-kit-stats kit)))
    (if (plist-get s :sex)
        (let ((sex (plist-get s :sex)))
          (list :x-het-rate (plist-get sex :x-het-rate)
                :y-called-rate (plist-get sex :y-call-rate)
                :method (plist-get sex :method)))
      (let* ((xhet (or (plist-get s :x-het) 0))
             (xhom (or (plist-get s :x-hom) 0))
             (xhemi (or (plist-get s :x-hemi) 0))
             (diploid (+ xhet xhom))
             (ytotal (or (plist-get s :y-total) 0)))
        (list :x-het-rate (cond ((> diploid 0) (/ (float xhet) diploid))
                                ((> xhemi 0) 0.0))
              :x-calls (+ diploid xhemi)
              :y-called (plist-get s :y-called)
              :y-called-rate (when (> ytotal 0)
                               (/ (float (plist-get s :y-called)) ytotal))
              :method (if (eq (genetics-kit-ref-calls kit) 'absent-means-ref)
                          "variant-only: non-PAR X heterozygosity among variant sites"
                        "non-PAR X heterozygosity + Y call rate"))))))

(defun genetics-infer-sex (kit)
  "Infer chromosomal gender of KIT: `male', `female' or `uncertain'.
X heterozygosity excludes the pseudoautosomal regions (see
`genetics-par-regions').  Arrays and other kits that list every assayed
site: male when at least half of the Y sites are called and X
heterozygosity is below 5%; female when X heterozygosity is at least 5%
and fewer than 20% of Y sites are called.  Variant-only WGS VCFs list
only called sites, so the Y call rate means nothing: male when X
heterozygosity among variant sites is below 15% and Y has calls, female
at 30% or more.  Anything else, and offset-indexed kits, is uncertain."
  (let ((cli (plist-get (plist-get (genetics-kit-stats kit) :sex) :call)))
    (cond
     (cli (intern cli))
     ((genetics-kit-lazy kit) 'uncertain)
     (t
      (let* ((e (genetics-sex-evidence kit))
             (xr (plist-get e :x-het-rate))
             (yr (plist-get e :y-called-rate)))
        (cond ((null xr) 'uncertain)
              ((eq (genetics-kit-ref-calls kit) 'absent-means-ref)
               (cond ((>= xr 0.30) 'female)
                     ((and (< xr 0.15) (> (or (plist-get e :y-called) 0) 0))
                      'male)
                     (t 'uncertain)))
              ((and (< xr 0.05) yr (>= yr 0.5)) 'male)
              ((and (>= xr 0.05) (or (null yr) (< yr 0.2))) 'female)
              (t 'uncertain)))))))

(defun genetics-sex-description (kit)
  "Return a sentence describing the gender inference for KIT."
  (let* ((e (genetics-sex-evidence kit))
         (xr (plist-get e :x-het-rate))
         (yr (plist-get e :y-called-rate)))
    (format "%s (non-PAR X heterozygosity %s, Y called %s)"
            (genetics-infer-sex kit)
            (if xr (format "%.2f%%" (* 100 xr)) "n/a")
            (if yr (format "%.0f%%" (* 100 yr)) "n/a"))))

(defun genetics-kit-rsid-free-p (kit)
  "Return non-nil if KIT has no rsids (records are keyed by position)."
  (let ((rsids (plist-get (genetics-kit-stats kit) :rsids)))
    (if rsids
        (zerop rsids)
      (eq (genetics-kit-has-rsids kit) 'none))))

(defun genetics-kit-caveats (kit)
  "Return a list of caveat strings for KIT (strand, build, chip, ref model)."
  (delq nil
        (append
         (list
          (genetics-kit-strand-note kit)
          (format "Positions are on GRCh%s; they cannot be compared to data on another build without liftover, which this package does not do."
                  (or (genetics-kit-build kit) "(unknown build)"))
          (when (equal (genetics-kit-chip kit) "unknown")
            (when (eq (genetics-kit-format kit) '23andme)
              "The 23andMe chip version could not be determined."))
          (when (eq (genetics-kit-assay kit) 'array)
            "A genotyping array reads a fixed subset of known sites (about 0.02% of the genome); a site that is not on the chip is simply absent.")
          (when (genetics-kit-rsid-free-p kit)
            (format "This file has no rsids; curated SNPs are matched by their GRCh%s position (coordinates verified against dbSNP and Ensembl, cited in the annotation file)."
                    (or (genetics-kit-build kit) "?")))
          (pcase (genetics-kit-ref-calls kit)
            ('absent-means-ref "Variant-only whole-genome VCF: a curated site that is absent is shown as homozygous reference, labelled inferred. It was not observed; an uncovered site would look the same.")
            ('unknown (when (eq (genetics-kit-format kit) 'vcf)
                        "This VCF may list only variant sites; absent sites are reported as not present, never as homozygous reference. Pass :ref-calls to `genetics-parse-file' or set `genetics-vcf-ref-calls' if the file is a whole-genome variant-only VCF.")))
          "Palindromic SNPs (A/T and C/G) cannot have their strand verified from the genotype alone.")
         (plist-get (genetics-kit-stats kit) :caveats)
         (genetics-kit-warnings kit))))

(defun genetics-builds-differ-p (a b)
  "Return non-nil if kits A and B have different known builds."
  (let ((ba (genetics-kit-build a)) (bb (genetics-kit-build b)))
    (and ba bb (not (equal ba bb)))))

(provide 'genetics-stats)
;;; genetics-stats.el ends here
