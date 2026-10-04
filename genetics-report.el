;;; genetics-report.el --- Org report generator for genetics.el -*- lexical-binding: t; -*-

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

;; `genetics-report' builds an Org buffer describing a kit: metadata,
;; strand and build caveats, a table of annotated hits, the APOE haplotype
;; and a disclaimer.  It can also be written to a file.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'org)
(require 'genetics-core)
(require 'genetics-stats)
(require 'genetics-annotate)

(defconst genetics-report-disclaimer
  "This report is informational only and is not medical advice. Consumer genotyping arrays and consumer sequencing are not diagnostic and can produce false positives; arrays cover only a small part of the genome. Discuss any finding with a qualified clinician or genetic counselor and confirm it with a clinical-grade test before acting on it."
  "Disclaimer text placed at the end of every report.")

(defun genetics--org-cell (value)
  "Return VALUE formatted as a safe Org table cell."
  (let ((s (if value (format "%s" value) "")))
    (replace-regexp-in-string
     "[ \t\n]+" " " (replace-regexp-in-string "|" "/" s))))

(defun genetics--report-hit-row (ann snp)
  "Return the Org table row string for annotation ANN and call SNP."
  (let* ((a (genetics-assess-snp ann snp))
         (url (genetics-annotation-url ann)))
    (format "| %s |\n"
            (mapconcat
             #'genetics--org-cell
             (list (genetics-annotation-rsid ann)
                   (genetics-annotation-gene ann)
                   (genetics-snp-genotype-label snp)
                   (genetics-annotation-risk-allele ann)
                   (let ((c (plist-get a :copies))) (if c c "n/a"))
                   (plist-get a :interpretation)
                   (genetics-annotation-strand ann)
                   (if url (format "[[%s][source]]" url) ""))
             " | "))))

(defun genetics-report-string (kit)
  "Return the Org report text for KIT."
  (let* ((stats (genetics-kit-stats kit))
         (hits (genetics-annotated-snps kit))
         (apoe (genetics-apoe-for-kit kit)))
    (with-temp-buffer
      (insert (format "#+TITLE: Genetics report: %s\n" (genetics-kit-name kit)))
      (insert (format "#+DATE: %s\n\n" (format-time-string "%Y-%m-%d")))
      (insert "* Kit\n")
      (insert (format "- File: %s\n" (abbreviate-file-name (genetics-kit-file kit))))
      (insert (format "- Format: %s\n" (genetics--format-label (genetics-kit-format kit))))
      (insert (format "- Build: GRCh%s\n" (or (genetics-kit-build kit) "unknown")))
      (insert (format "- Chip: %s\n" (or (genetics-kit-chip kit) "unknown")))
      (when (genetics-kit-sample kit)
        (insert (format "- Sample: %s\n" (genetics-kit-sample kit))))
      (insert (format "- Records: %d\n" (plist-get stats :total)))
      (when (genetics-kit-assay kit)
        (insert (format "- Assay: %s (reference calls: %s)\n"
                        (genetics-kit-assay kit)
                        (or (genetics-kit-ref-calls kit) "unknown"))))
      (when (plist-get stats :nocalls)
        (insert (format "- No-calls: %d (%.2f%%)\n" (plist-get stats :nocalls)
                        (* 100 (genetics-nocall-rate kit)))))
      (unless (and (genetics-kit-lazy kit) (not (plist-get stats :sex)))
        (insert (format "- Inferred sex: %s\n" (genetics-sex-description kit))))
      (insert "\n* Strand and build caveats\n")
      (dolist (c (genetics-kit-caveats kit))
        (insert (format "- %s\n" c)))
      (insert "- Risk alleles in this report are given on the + strand (identical on GRCh37 and GRCh38 for the curated SNPs). A genotype containing only complementary alleles is flagged as a possible strand flip and is never flipped automatically.\n")
      (insert "\n* Annotated findings\n")
      (if (null hits)
          (insert "No annotated SNPs were found in this kit.\n")
        (insert "| rsid | gene | genotype | risk allele | copies | interpretation | strand note | source |\n|-\n")
        (dolist (h hits)
          (insert (genetics--report-hit-row (car h) (cdr h))))
        (when (cl-some (lambda (h) (genetics-snp-inferred-p (cdr h))) hits)
          (insert "\nGenotypes marked \"(inferred ref)\" were not observed. " genetics-inferred-ref-text "\n")))
      (insert "\n* APOE haplotype\n")
      (if (null apoe)
          (insert "rs429358 and rs7412 were not both present, so APOE could not be assessed.\n")
        (insert (format "- Diplotype: %s\n- Status: %s\n- %s\n"
                        (or (plist-get apoe :diplotype) "undetermined")
                        (plist-get apoe :status)
                        (plist-get apoe :description)))
        (insert "- Method: rs429358 T/C and rs7412 C/T on the + strand; e2 = T+T, e3 = T+C, e4 = C+C.\n"))
      (insert "\n* Disclaimer\n" genetics-report-disclaimer "\n")
      (buffer-string))))

(defun genetics-report-write (kit file)
  "Write the Org report for KIT to FILE and return FILE."
  (let ((text (genetics-report-string kit)))
    (with-temp-file file (insert text))
    file))

;;;###autoload
(defun genetics-report (&optional kit file)
  "Show an Org report for KIT; with FILE (or a prefix argument) also save it.
Returns the report buffer."
  (interactive
   (let ((kit (genetics--read-kit "Report for kit: ")))
     (list kit (when current-prefix-arg
                 (read-file-name "Write report to: ")))))
  (let* ((kit (genetics-find-kit (or kit (genetics--read-kit "Report for kit: "))))
         (buf (get-buffer-create (format "*genetics-report: %s*"
                                         (genetics-kit-name kit))))
         (text (genetics-report-string kit)))
    (with-current-buffer buf
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert text)
        (org-mode)
        (goto-char (point-min))))
    (when file
      (with-temp-file file (insert text)))
    (pop-to-buffer buf)
    buf))

(provide 'genetics-report)
;;; genetics-report.el ends here
