;;; genetics-export.el --- Export genetics views to CSV and JSON -*- lexical-binding: t; -*-

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

;; Writes the browser's filtered view (or a whole kit) to CSV or JSON.  The
;; CSV header is RSID,CHROMOSOME,POSITION,RESULT plus extra columns, so an
;; exported file can be parsed again by `genetics-open'.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)
(require 'genetics-core)
(require 'genetics-annotate)
(require 'genetics-browse)

(defun genetics-export--context (kit filters)
  "Return (KIT . FILTERS) from KIT, FILTERS or the current browser buffer."
  (if (and (null kit) (derived-mode-p 'genetics-browse-mode))
      (cons genetics--buffer-kit genetics-browse--filters)
    (cons (genetics-find-kit (or kit (genetics--read-kit "Export kit: ")))
          filters)))

(defun genetics-export--snps (kit filters)
  "Return all SNPs of KIT matching FILTERS (no display limit)."
  (car (genetics-browse-rows kit filters nil)))

(defun genetics-export--csv-quote (value)
  "Return VALUE as a quoted CSV field."
  (format "\"%s\"" (replace-regexp-in-string "\"" "\"\"" (format "%s" (or value "")))))

;;;###autoload
(defun genetics-export-csv (file &optional kit filters)
  "Write records to CSV FILE and return the number written.
In a browser buffer the buffer's kit and filters are used; otherwise KIT
and plist FILTERS (see `genetics-browse-rows') select the records."
  (interactive (list (read-file-name "Export CSV to: ")))
  (let* ((ctx (genetics-export--context kit filters))
         (snps (genetics-export--snps (car ctx) (cdr ctx)))
         (annotations (genetics-annotations)))
    (genetics--with-output-file file
      (insert "RSID,CHROMOSOME,POSITION,RESULT,ZYGOSITY,GENE,REF,ALT\n")
      (dolist (s snps)
        (let ((ann (genetics-snp-annotation s (genetics-kit-build (car ctx))
                                            annotations)))
          (insert (mapconcat
                   #'genetics-export--csv-quote
                   (list (genetics-snp-rsid s) (genetics-snp-chrom s)
                         (genetics-snp-pos s) (genetics-snp-genotype s)
                         (genetics-zygosity (genetics-snp-genotype s))
                         (and ann (genetics-annotation-gene ann))
                         (genetics-snp-ref s) (genetics-snp-alt s))
                   ",")
                  "\n"))))
    (length snps)))

;;;###autoload
(defun genetics-export-json (file &optional kit filters)
  "Write records to JSON FILE and return the number written.
Arguments FILE, KIT and FILTERS are as for `genetics-export-csv'."
  (interactive (list (read-file-name "Export JSON to: ")))
  (let* ((ctx (genetics-export--context kit filters))
         (snps (genetics-export--snps (car ctx) (cdr ctx)))
         (annotations (genetics-annotations)))
    (genetics--with-output-file file
      (insert
       (genetics--json-serialize
        (vconcat
         (mapcar
          (lambda (s)
            (let ((ann (genetics-snp-annotation
                        s (genetics-kit-build (car ctx)) annotations)))
              `((rsid . ,(genetics-snp-rsid s))
                (chromosome . ,(genetics-snp-chrom s))
                (position . ,(genetics-snp-pos s))
                (genotype . ,(genetics-snp-genotype s))
                (zygosity . ,(symbol-name
                              (genetics-zygosity (genetics-snp-genotype s))))
                (gene . ,(or (and ann (genetics-annotation-gene ann)) :null))
                (ref . ,(or (genetics-snp-ref s) :null))
                (alt . ,(or (genetics-snp-alt s) :null)))))
          snps)))
       "\n"))
    (length snps)))

(provide 'genetics-export)
;;; genetics-export.el ends here
