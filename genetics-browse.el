;;; genetics-browse.el --- Tabulated browser for genetics.el kits -*- lexical-binding: t; -*-

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

;; `genetics-browse' lists the records of a kit in a tabulated buffer with
;; filters (chromosome, position range, rsid regexp, genotype, zygosity,
;; annotated only).  Display is capped by `genetics-browse-limit'.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'tabulated-list)
(require 'genetics-core)
(require 'genetics-annotate)

(declare-function genetics-lookup "genetics-lookup" (rsid))
(declare-function genetics-report "genetics-report" (&optional kit file))
(declare-function genetics-export-csv "genetics-export" (file &optional kit filters))
(declare-function genetics-export-json "genetics-export" (file &optional kit filters))
(declare-function genetics-summary "genetics" (&optional kit))

(defvar-local genetics-browse--filters nil
  "Plist of active filters in the current browser buffer.
Keys: :chrom :range (START . END) :rsid :genotype :zygosity :annotated.")

(defvar-local genetics-browse--truncated nil
  "Non-nil when the last refresh hid rows because of the display limit.")

(defun genetics-browse-snp-matches-p (snp filters annotations &optional build)
  "Return non-nil if SNP satisfies plist FILTERS.
ANNOTATIONS is the annotation hash table; BUILD, the kit's build, lets
records without an rsid match an annotation by position."
  (let ((chrom (plist-get filters :chrom))
        (range (plist-get filters :range))
        (rx (plist-get filters :rsid))
        (gt (plist-get filters :genotype))
        (zyg (plist-get filters :zygosity)))
    (and (or (null chrom) (equal chrom (genetics-snp-chrom snp)))
         (or (null range)
             (<= (car range) (genetics-snp-pos snp) (cdr range)))
         (or (null rx)
             (let ((case-fold-search t))
               (string-match-p rx (genetics-snp-rsid snp))))
         (or (null gt)
             (equal gt (genetics-genotype-key (genetics-snp-genotype snp))))
         (or (null zyg)
             (eq zyg (genetics-zygosity (genetics-snp-genotype snp))))
         (or (null (plist-get filters :annotated))
             (genetics-snp-annotation snp build annotations)))))

(defun genetics-browse-rows (kit filters &optional limit)
  "Return (SNPS . TRUNCATED) for KIT filtered by FILTERS.
At most LIMIT records are returned (nil means all); TRUNCATED is non-nil
when more records matched than were returned."
  (let ((annotations (genetics-annotations)) (rows nil) (n 0) (truncated nil))
    (genetics-kit-map-snps
     kit
     (lambda (snp)
       (when (genetics-browse-snp-matches-p snp filters annotations
                                            (genetics-kit-build kit))
         (if (and limit (>= n limit))
             (progn (setq truncated t) 'stop)
           (push snp rows)
           (cl-incf n)
           nil)))
     (plist-get filters :chrom))
    (cons (nreverse rows) truncated)))

(defun genetics-browse-describe-filters (filters)
  "Return a short string describing FILTERS."
  (let ((parts nil))
    (when (plist-get filters :chrom)
      (push (format "chr=%s" (plist-get filters :chrom)) parts))
    (when (plist-get filters :range)
      (push (format "pos=%d-%d" (car (plist-get filters :range))
                    (cdr (plist-get filters :range)))
            parts))
    (when (plist-get filters :rsid)
      (push (format "rsid~%s" (plist-get filters :rsid)) parts))
    (when (plist-get filters :genotype)
      (push (format "genotype=%s" (plist-get filters :genotype)) parts))
    (when (plist-get filters :zygosity)
      (push (format "%s" (plist-get filters :zygosity)) parts))
    (when (plist-get filters :annotated) (push "annotated" parts))
    (if parts (string-join (nreverse parts) " ") "none")))

(defun genetics-browse--entry (snp annotations &optional build)
  "Return a tabulated-list entry for SNP using ANNOTATIONS on BUILD."
  (let ((ann (genetics-snp-annotation snp build annotations)))
    (list (genetics-snp-rsid snp)
          (vector (genetics-snp-rsid snp) (genetics-snp-chrom snp)
                  (number-to-string (genetics-snp-pos snp))
                  (genetics-snp-genotype snp)
                  (symbol-name (genetics-zygosity (genetics-snp-genotype snp)))
                  (or (and ann (genetics-annotation-gene ann)) "")))))

(defun genetics-browse--entries ()
  "Return the entries for the current buffer, updating the header line."
  (let* ((res (genetics-browse-rows genetics--buffer-kit genetics-browse--filters
                                    genetics-browse-limit))
         (annotations (genetics-annotations)))
    (setq genetics-browse--truncated (cdr res))
    (setq header-line-format
          (format " Filters: %s | %d rows%s"
                  (genetics-browse-describe-filters genetics-browse--filters)
                  (length (car res))
                  (if (cdr res)
                      (format " (truncated at `genetics-browse-limit' = %d)"
                              genetics-browse-limit)
                    "")))
    (when (cdr res)
      (message "Showing first %d matches; refine filters or raise `genetics-browse-limit'"
               genetics-browse-limit))
    (mapcar (lambda (s) (genetics-browse--entry
                         s annotations (genetics-kit-build genetics--buffer-kit)))
            (car res))))

(defun genetics-browse--pos< (a b)
  "Return non-nil if entry A has a smaller position than entry B."
  (< (string-to-number (aref (cadr a) 2)) (string-to-number (aref (cadr b) 2))))

(defvar genetics-browse-mode-map
  (let ((map (make-sparse-keymap)))
    (set-keymap-parent map tabulated-list-mode-map)
    (define-key map (kbd "RET") #'genetics-browse-lookup)
    (define-key map "c" #'genetics-browse-filter-chromosome)
    (define-key map "p" #'genetics-browse-filter-range)
    (define-key map "s" #'genetics-browse-filter-rsid)
    (define-key map "g" #'genetics-browse-filter-genotype)
    (define-key map "n" #'genetics-browse-filter-no-calls)
    (define-key map "h" #'genetics-browse-filter-heterozygous)
    (define-key map "o" #'genetics-browse-filter-homozygous)
    (define-key map "a" #'genetics-browse-filter-annotated)
    (define-key map "x" #'genetics-browse-clear-filters)
    (define-key map "E" #'genetics-browse-export-csv)
    (define-key map "J" #'genetics-browse-export-json)
    (define-key map "R" #'genetics-browse-report)
    (define-key map "i" #'genetics-browse-summary)
    map)
  "Keymap for `genetics-browse-mode'.")

(define-derived-mode genetics-browse-mode tabulated-list-mode "Genetics-Browse"
  "Major mode listing the records of a genetics kit.

\\{genetics-browse-mode-map}"
  (setq tabulated-list-format
        [("rsid" 14 t) ("chrom" 6 t) ("position" 11 genetics-browse--pos<)
         ("genotype" 10 t) ("zygosity" 13 t) ("gene" 14 t)])
  (setq tabulated-list-padding 1)
  ;; The header line shows the filters, so the column headings go in the
  ;; buffer's first line instead.
  (setq tabulated-list-use-header-line nil)
  (setq tabulated-list-entries #'genetics-browse--entries)
  (tabulated-list-init-header))

(defun genetics-browse--require-buffer ()
  "Signal a user error unless in a browser buffer."
  (unless (derived-mode-p 'genetics-browse-mode)
    (user-error "Not in a genetics browser buffer")))

(defun genetics-browse--set (key value)
  "Set filter KEY to VALUE (nil clears it) and refresh the buffer."
  (genetics-browse--require-buffer)
  (setq genetics-browse--filters (plist-put genetics-browse--filters key value))
  (tabulated-list-revert)
  genetics-browse--filters)

;;;###autoload
(defun genetics-browse (&optional kit)
  "Browse the records of KIT (default: ask) in a filterable table."
  (interactive (list (genetics--read-kit "Browse kit: ")))
  (let* ((kit (genetics-find-kit (or kit (genetics--read-kit "Browse kit: "))))
         (buf (get-buffer-create (format "*genetics-browse: %s*"
                                         (genetics-kit-name kit)))))
    (with-current-buffer buf
      (genetics-browse-mode)
      (setq genetics--buffer-kit kit)
      (setq genetics-browse--filters nil)
      (tabulated-list-revert))
    (pop-to-buffer buf)
    buf))

(defun genetics-browse-filter-chromosome (chrom)
  "Show only chromosome CHROM (empty string clears the filter)."
  (interactive (list (read-string "Chromosome (empty to clear): ")))
  (genetics-browse--set
   :chrom (unless (string-empty-p chrom) (genetics-normalize-chrom chrom))))

(defun genetics-browse-filter-range (start end)
  "Show only positions from START to END (inclusive)."
  (interactive (list (read-number "Start position: ")
                     (read-number "End position: ")))
  (when (> start end)
    (user-error "Start position %d is after end %d" start end))
  (genetics-browse--set :range (cons start end)))

(defun genetics-browse-filter-rsid (regexp)
  "Show only records whose rsid matches REGEXP (empty clears)."
  (interactive (list (read-regexp "rsid regexp (empty to clear): " "")))
  (unless (string-empty-p regexp) (string-match-p regexp ""))
  (genetics-browse--set :rsid (unless (string-empty-p regexp) regexp)))

(defun genetics-browse-filter-genotype (genotype)
  "Show only records with GENOTYPE, ignoring allele order (empty clears)."
  (interactive (list (read-string "Genotype (empty to clear): ")))
  (genetics-browse--set
   :genotype (unless (string-empty-p genotype)
               (genetics-genotype-key (upcase genotype)))))

(defun genetics-browse-filter-no-calls ()
  "Show only no-calls."
  (interactive)
  (genetics-browse--set :zygosity 'no-call))

(defun genetics-browse-filter-heterozygous ()
  "Show only heterozygous records."
  (interactive)
  (genetics-browse--set :zygosity 'heterozygous))

(defun genetics-browse-filter-homozygous ()
  "Show only homozygous records."
  (interactive)
  (genetics-browse--set :zygosity 'homozygous))

(defun genetics-browse-filter-annotated ()
  "Toggle showing only annotated records."
  (interactive)
  (genetics-browse--set :annotated (not (plist-get genetics-browse--filters
                                                   :annotated))))

(defun genetics-browse-clear-filters ()
  "Remove all filters."
  (interactive)
  (genetics-browse--require-buffer)
  (setq genetics-browse--filters nil)
  (tabulated-list-revert))

(defun genetics-browse-lookup ()
  "Show the detail buffer for the record at point."
  (interactive)
  (genetics-browse--require-buffer)
  (let ((id (tabulated-list-get-id)))
    (unless id (user-error "No record on this line"))
    (require 'genetics-lookup)
    (genetics-lookup id)))

(defun genetics-browse-export-csv (file)
  "Export the filtered view to CSV FILE."
  (interactive (list (read-file-name "Export CSV to: ")))
  (require 'genetics-export)
  (genetics-export-csv file))

(defun genetics-browse-export-json (file)
  "Export the filtered view to JSON FILE."
  (interactive (list (read-file-name "Export JSON to: ")))
  (require 'genetics-export)
  (genetics-export-json file))

(defun genetics-browse-report ()
  "Open the report for the browsed kit."
  (interactive)
  (require 'genetics-report)
  (genetics-report genetics--buffer-kit))

(defun genetics-browse-summary ()
  "Open the summary for the browsed kit."
  (interactive)
  (require 'genetics)
  (genetics-summary genetics--buffer-kit))

(provide 'genetics-browse)
;;; genetics-browse.el ends here
