;;; genetics-lookup.el --- Single-SNP lookup buffer for genetics.el -*- lexical-binding: t; -*-

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

;; `genetics-lookup' shows the genotype of one rsid in every loaded kit,
;; with its annotation and interpretation, and (opt-in) SNPedia text.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'button)
(require 'genetics-core)
(require 'genetics-annotate)

(declare-function genetics-snpedia-summary "genetics-snpedia" (rsid))

(defvar genetics-lookup-mode-map
  (let ((map (make-sparse-keymap)))
    (set-keymap-parent map special-mode-map)
    map)
  "Keymap for `genetics-lookup-mode'.")

(define-derived-mode genetics-lookup-mode special-mode "Genetics-Lookup"
  "Major mode for the single-SNP detail buffer.

\\{genetics-lookup-mode-map}")

(defun genetics--normalize-rsid (rsid)
  "Return RSID trimmed and lower-cased."
  (downcase (string-trim rsid)))

(defun genetics--all-ids ()
  "Return a hash table of ids across all loaded kits."
  (if (null (cdr genetics-loaded-kits))
      (genetics-kit-ids (car genetics-loaded-kits))
    (let ((all (make-hash-table :test 'equal)))
      (dolist (k genetics-loaded-kits)
        (maphash (lambda (id _) (puthash id t all)) (genetics-kit-ids k)))
      all)))

(defun genetics--insert-url (url)
  "Insert URL as a browse button."
  (insert-text-button url 'action (lambda (b) (browse-url (button-label b)))
                      'follow-link t))

(defun genetics--lookup-annotation-section (ann calls)
  "Insert the annotation section for ANN.
CALLS is a list of (KIT-NAME . SNP) for the kits that have the record."
  (insert (propertize "Annotation\n" 'face 'bold))
  (cl-flet ((row (label value)
              (when value (insert (format "  %-12s %s\n" label value)))))
    (row "Gene:" (genetics-annotation-gene ann))
    (row "Risk allele:" (genetics-annotation-risk-allele ann))
    (row "Magnitude:" (genetics-annotation-magnitude ann))
    (row "Effect:" (genetics-annotation-effect ann))
    (row "Strand:" (genetics-annotation-strand ann))
    (row "Notes:" (genetics-annotation-notes ann))
    (dolist (call calls)
      (let ((a (genetics-assess-snp ann (cdr call))))
        (insert (format "\n  %s: %s\n" (car call)
                        (genetics-snp-genotype-label (cdr call))))
        (row "Copies:" (plist-get a :copies))
        (row "Flag:" (plist-get a :flag))
        (row "Meaning:" (plist-get a :interpretation)))))
  (when (genetics-annotation-url ann)
    (insert "\n  Source:      ")
    (genetics--insert-url (genetics-annotation-url ann))
    (insert "\n")))

(defun genetics--lookup-apoe-section (rsid)
  "Insert APOE haplotype text when RSID is one of the APOE SNPs."
  (when (member rsid genetics-apoe-snps)
    (dolist (kit genetics-loaded-kits)
      (let ((res (genetics-apoe-for-kit kit)))
        (when res
          (insert (format "\nAPOE haplotype in %s: %s\n  %s\n"
                          (genetics-kit-name kit)
                          (or (plist-get res :diplotype) "undetermined")
                          (plist-get res :description))))))))

;;;###autoload
(defun genetics-lookup (rsid)
  "Show genotype, annotation and interpretation of RSID in loaded kits.
Returns the lookup buffer."
  (interactive
   (progn
     (unless genetics-loaded-kits
       (genetics--error 'genetics-no-kit
                        "No kit loaded; run M-x genetics-open first"))
     (list (completing-read "rsid: " (genetics--all-ids)))))
  (unless genetics-loaded-kits
    (genetics--error 'genetics-no-kit
                     "No kit loaded; run M-x genetics-open first"))
  (let* ((rsid (genetics--normalize-rsid rsid))
         (ann (genetics-annotation rsid))
         (buf (get-buffer-create (format "*genetics-lookup: %s*" rsid)))
         (calls nil))
    (with-current-buffer buf
      (genetics-lookup-mode)
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert (propertize (format "%s\n\n" rsid) 'face 'bold))
        (dolist (kit genetics-loaded-kits)
          (let ((snp (genetics-kit-resolve kit rsid)))
            (when snp (push (cons (genetics-kit-name kit) snp) calls))
            (insert
             (if snp
                 (format "%-24s %s  chr%s:%d  %s (GRCh%s)%s\n"
                         (genetics-kit-name kit)
                         (genetics-snp-genotype-label snp)
                         (genetics-snp-chrom snp) (genetics-snp-pos snp)
                         (genetics-zygosity (genetics-snp-genotype snp))
                         (or (genetics-kit-build kit) "?")
                         (cond ((genetics-snp-inferred-p snp)
                                "\n                         inferred, not observed: absent from a variant-only WGS VCF")
                               ((not (equal (genetics-snp-rsid snp) rsid))
                                "  matched by position")
                               (t "")))
               (format "%-24s not present\n" (genetics-kit-name kit))))))
        (insert "\n")
        ;; a position id (chrom:pos) from a kit without rsids
        (unless ann
          (cl-loop for kit in genetics-loaded-kits
                   for snp = (cdr (assoc (genetics-kit-name kit) calls))
                   until ann
                   when snp
                   do (setq ann (genetics-snp-annotation
                                 snp (genetics-kit-build kit)))))
        (if ann
            (genetics--lookup-annotation-section ann (nreverse calls))
          (insert "No annotation for this rsid in `genetics-annotation-files'.\n"))
        (genetics--lookup-apoe-section
         (if ann (genetics-annotation-rsid ann) rsid))
        (when genetics-snpedia-enabled
          (require 'genetics-snpedia)
          (insert "\nSNPedia\n")
          (condition-case err
              (let ((s (genetics-snpedia-summary rsid)))
                (insert (format "  Gene: %s\n  Orientation: %s\n  Summary: %s\n"
                                (or (plist-get s :gene) "n/a")
                                (or (plist-get s :orientation) "n/a")
                                (or (plist-get s :summary) "n/a"))))
            (genetics-error
             (insert (format "  SNPedia lookup failed: %s\n"
                             (error-message-string err))))))
        (goto-char (point-min))))
    (pop-to-buffer buf)
    buf))

(provide 'genetics-lookup)
;;; genetics-lookup.el ends here
