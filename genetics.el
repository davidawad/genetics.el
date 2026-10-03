;;; genetics.el --- Read and explore consumer genetics raw-data files -*- lexical-binding: t; -*-

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

;; Open 23andMe, AncestryDNA, MyHeritage/FTDNA CSV and VCF files, browse
;; and filter genotypes, annotate them from editable JSON/Org files,
;; generate Org reports and compare kits -- all locally.
;;
;;   M-x genetics-open      parse a file, show its summary
;;   M-x genetics-browse    filterable table of records
;;   M-x genetics-lookup    one rsid across all loaded kits
;;   M-x genetics-report    Org report of annotated findings and APOE
;;   M-x genetics-compare   concordance between two kits
;;
;; PRIVACY: no data leaves your machine.  The only network code lives in
;; genetics-snpedia.el, is off by default, and sends only an rsid.
;; Results are informational, not medical advice.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'button)
(require 'genetics-core)
(require 'genetics-parse)
(require 'genetics-stats)
(require 'genetics-annotate)
(require 'genetics-browse)
(require 'genetics-lookup)
(require 'genetics-report)
(require 'genetics-compare)
(require 'genetics-export)

(defvar genetics-summary-mode-map
  (let ((map (make-sparse-keymap)))
    (set-keymap-parent map special-mode-map)
    (define-key map "b" #'genetics-summary-browse)
    (define-key map "r" #'genetics-summary-report)
    (define-key map "l" #'genetics-lookup)
    (define-key map "c" #'genetics-compare)
    (define-key map "g" #'genetics-summary-refresh)
    map)
  "Keymap for `genetics-summary-mode'.")

(define-derived-mode genetics-summary-mode special-mode "Genetics-Summary"
  "Major mode for the kit summary buffer.

\\{genetics-summary-mode-map}")

(defun genetics--insert-button (label action)
  "Insert a button LABEL running ACTION (a function of no arguments)."
  (insert-text-button label 'action (lambda (_b) (funcall action))
                      'follow-link t))

(defun genetics-summary-string (kit)
  "Return the plain-text summary of KIT."
  (let* ((s (genetics-kit-stats kit))
         (total (plist-get s :total)))
    (with-temp-buffer
      (insert (format "Kit:       %s\n" (genetics-kit-name kit)))
      (insert (format "File:      %s\n" (abbreviate-file-name (genetics-kit-file kit))))
      (insert (format "Format:    %s\n" (genetics--format-label (genetics-kit-format kit))))
      (insert (format "Build:     GRCh%s\n" (or (genetics-kit-build kit) "unknown")))
      (insert (format "Chip:      %s\n" (or (genetics-kit-chip kit) "unknown")))
      (when (genetics-kit-sample kit)
        (insert (format "Sample:    %s\n" (genetics-kit-sample kit))))
      (insert (format "SNPs:      %d%s\n" total
                      (if (genetics-kit-lazy kit)
                          " (offset-indexed; records are read on demand)" "")))
      (if (genetics-kit-lazy kit)
          (insert "No-calls, het/hom counts and sex inference are not computed in offset-indexed mode.\n")
        (insert (format "No-calls:  %d (%.2f%%)\n" (plist-get s :nocalls)
                        (* 100 (or (genetics-nocall-rate kit) 0))))
        (insert (format "Het/Hom:   %d heterozygous, %d homozygous, %d hemizygous\n"
                        (plist-get s :het) (plist-get s :hom)
                        (plist-get s :hemi)))
        (insert (format "Sex:       %s\n" (genetics-sex-description kit))))
      (insert "\nPer-chromosome counts\n")
      (dolist (c (plist-get s :chrom-counts))
        (insert (format "  %-4s %8d\n" (car c) (cdr c))))
      (insert "\nCaveats\n")
      (dolist (c (genetics-kit-caveats kit))
        (insert (format "  - %s\n" c)))
      (insert "  - Informational only; not medical advice.\n")
      (buffer-string))))

;;;###autoload
(defun genetics-summary (&optional kit)
  "Show the summary buffer for KIT (default: current or ask)."
  (interactive)
  (let* ((kit (genetics-find-kit (or kit (genetics--read-kit "Summary for kit: "))))
         (buf (get-buffer-create (format "*genetics: %s*" (genetics-kit-name kit)))))
    (with-current-buffer buf
      (genetics-summary-mode)
      (setq genetics--buffer-kit kit)
      (genetics-summary-refresh))
    (pop-to-buffer buf)
    buf))

(defun genetics-summary-refresh ()
  "Redraw the summary buffer."
  (interactive)
  (let ((kit genetics--buffer-kit) (inhibit-read-only t))
    (erase-buffer)
    (insert (genetics-summary-string kit) "\n")
    (genetics--insert-button "[Browse]" (lambda () (genetics-browse kit)))
    (insert " ")
    (genetics--insert-button "[Report]" (lambda () (genetics-report kit)))
    (insert " ")
    (genetics--insert-button "[Lookup]" #'genetics-lookup-interactively)
    (insert "\n")
    (goto-char (point-min))))

(defun genetics-lookup-interactively ()
  "Prompt for an rsid and look it up."
  (call-interactively #'genetics-lookup))

(defun genetics-summary-browse ()
  "Browse the kit of this summary buffer."
  (interactive)
  (genetics-browse genetics--buffer-kit))

(defun genetics-summary-report ()
  "Open the report of this summary buffer's kit."
  (interactive)
  (genetics-report genetics--buffer-kit))

;;;###autoload
(defun genetics-open (file)
  "Parse genotype FILE, register the kit and show its summary.
Returns the kit."
  (interactive
   (list (read-file-name "Genetics file: " genetics-data-directory nil t)))
  (let ((kit (genetics-parse-file file)))
    (setf (genetics-kit-name kit) (genetics--unique-name
                                   (genetics-kit-name kit)
                                   (genetics-kit-file kit)))
    (genetics-register-kit kit)
    (genetics-summary kit)
    (when (cl-some (lambda (k) (and (not (eq k kit))
                                    (genetics-builds-differ-p k kit)))
                   genetics-loaded-kits)
      (display-warning 'genetics
                       "Loaded kits use different genome builds; do not compare positions across them (no liftover)."))
    kit))

;;;###autoload
(defun genetics-close (kit)
  "Unload KIT (a kit or its name) and kill its buffers."
  (interactive (list (genetics--read-kit "Close kit: ")))
  (let ((kit (genetics-find-kit kit)))
    (setq genetics-loaded-kits (delq kit genetics-loaded-kits))
    (dolist (b (buffer-list))
      (when (eq (buffer-local-value 'genetics--buffer-kit b) kit)
        (kill-buffer b)))))

(provide 'genetics)
;;; genetics.el ends here
