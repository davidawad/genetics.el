;;; genetics-catalog.el --- Annotation catalog files for genetics.el -*- lexical-binding: t; -*-

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

;; Loads the user-editable JSON and Org annotation files (see
;; genetics-annotate.el for their format), validates their coordinates,
;; caches them and looks annotations up by rsid or position.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)
(require 'genetics-core)

(cl-defstruct (genetics-annotation (:constructor genetics-annotation-create)
                                   (:copier nil))
  "Annotation of one SNP."
  rsid gene risk-allele other-allele effect magnitude notes url strand
  genotypes coordinates sources)

(defun genetics-build-key (build)
  "Return \"37\" or \"38\" for BUILD (e.g. \"GRCh38\", \"hg19\", 37), else nil."
  (let ((b (downcase (format "%s" (or build "")))))
    (cond ((string-match-p "\\`\\(grch\\)?37\\(\\.[0-9p]+\\)?\\'\\|\\`hg19\\'\\|\\`b37\\'" b) "37")
          ((string-match-p "\\`\\(grch\\)?38\\(\\.[0-9p]+\\)?\\'\\|\\`hg38\\'\\|\\`b38\\'" b) "38"))))

(defun genetics--ann-site (file rsid build chrom pos ref)
  "Validate one coordinate of RSID in FILE; return (BUILD-KEY . PLIST).
BUILD, CHROM, POS and REF are the raw values; POS may be a string."
  (let ((key (genetics-build-key build))
        (pos (if (stringp pos) (and (string-match-p "\\`[0-9]+\\'" pos)
                                    (string-to-number pos))
               pos)))
    (unless key
      (genetics--ann-error file "%s: unknown build %S (use GRCh37 or GRCh38)"
                           rsid build))
    (unless (and (stringp chrom) (not (string-empty-p chrom)))
      (genetics--ann-error file "%s: %s coordinate has no chrom" rsid build))
    (unless (and (integerp pos) (> pos 0))
      (genetics--ann-error file "%s: %s coordinate needs a positive pos" rsid
                           build))
    (when (and ref (not (string-match-p "\\`[ACGTacgt]+\\'" ref)))
      (genetics--ann-error file "%s: %s ref %S is not a base" rsid build ref))
    (cons key (list :chrom (genetics-normalize-chrom chrom) :pos pos
                    :ref (and ref (upcase ref))))))

;;;; Loading

(defun genetics--ann-error (file fmt &rest args)
  "Signal an annotation error for FILE with message FMT and ARGS."
  (genetics--error 'genetics-annotation-error "%s: %s" file
                   (apply #'format fmt args)))

(defun genetics--ann-genotypes (alist)
  "Convert ALIST of (GENOTYPE . TEXT) to normalized-key pairs."
  (delq nil
        (mapcar (lambda (c)
                  (let ((k (genetics-genotype-key
                            (upcase (format "%s" (car c))))))
                    (when k (cons k (cdr c)))))
                alist)))

(defun genetics--ann-from-json (obj file index)
  "Build an annotation from decoded JSON OBJ (item INDEX of FILE)."
  (let ((rsid (alist-get 'rsid obj)))
    (unless (and (stringp rsid) (not (string-empty-p rsid)))
      (genetics--ann-error file "entry %d has no rsid" index))
    (genetics-annotation-create
     :rsid (downcase rsid) :gene (alist-get 'gene obj)
     :risk-allele (let ((r (alist-get 'risk_allele obj))) (and r (upcase r)))
     :other-allele (let ((r (alist-get 'other_allele obj))) (and r (upcase r)))
     :effect (alist-get 'effect obj)
     :magnitude (alist-get 'magnitude obj)
     :notes (alist-get 'notes obj) :url (alist-get 'url obj)
     :strand (alist-get 'strand obj)
     :genotypes (genetics--ann-genotypes (alist-get 'genotypes obj))
     :coordinates
     (mapcar (lambda (c)
               (let ((v (cdr c)))
                 (unless (and (consp v) (consp (car v)))
                   (genetics--ann-error file "%s: coordinates.%s must be an object"
                                        rsid (car c)))
                 (genetics--ann-site file rsid (symbol-name (car c))
                                     (let ((ch (alist-get 'chrom v)))
                                       (if (numberp ch) (number-to-string ch) ch))
                                     (alist-get 'pos v) (alist-get 'ref v))))
             (alist-get 'coordinates obj))
     :sources (alist-get 'coordinate_sources obj))))

(defun genetics--load-annotation-json (file)
  "Return the annotations in JSON FILE."
  (let ((data (condition-case err
                  (with-temp-buffer
                    (let ((coding-system-for-read 'utf-8))
                      (insert-file-contents file))
                    (genetics--json-parse (buffer-string)
                                          :object-type 'alist :array-type 'list
                                          :null-object nil :false-object nil))
                (json-error
                 (genetics--ann-error file "invalid JSON (%s)"
                                      (error-message-string err))))))
    (unless (and (listp data)
                 (cl-every (lambda (o) (and (consp o) (consp (car o)))) data))
      (genetics--ann-error file "top level must be an array of objects"))
    (let ((i 0))
      (mapcar (lambda (o) (genetics--ann-from-json o file (cl-incf i))) data))))

(defun genetics--ann-from-org (props body file)
  "Build an annotation from Org PROPS (alist) and BODY text of FILE."
  (let ((rsid (cdr (assoc "RSID" props)))
        (get (lambda (k) (cdr (assoc k props)))))
    (unless rsid
      (genetics--ann-error file "heading without an RSID or rsNNN title"))
    (genetics-annotation-create
     :rsid (downcase rsid) :gene (funcall get "GENE")
     :risk-allele (let ((r (funcall get "RISK_ALLELE"))) (and r (upcase r)))
     :other-allele (let ((r (funcall get "OTHER_ALLELE"))) (and r (upcase r)))
     :effect (funcall get "EFFECT")
     :magnitude (let ((m (funcall get "MAGNITUDE")))
                  (and m (if (string-match-p "\\`[0-9.]+\\'" m)
                             (string-to-number m) m)))
     :notes (let ((n (or (funcall get "NOTES") (string-trim body))))
              (and (not (string-empty-p n)) n))
     :url (funcall get "URL") :strand (funcall get "STRAND")
     :genotypes (genetics--ann-genotypes
                 (cl-loop for (k . v) in props
                          when (string-prefix-p "GT_" k)
                          collect (cons (substring k 3) v)))
     :coordinates
     (cl-loop for b in '("GRCH37" "GRCH38")
              for chrom = (funcall get (concat b "_CHROM"))
              for pos = (funcall get (concat b "_POS"))
              when (or chrom pos)
              collect (genetics--ann-site file rsid b chrom pos
                                          (funcall get (concat b "_REF"))))
     :sources (let ((src (funcall get "COORDINATE_SOURCES")))
                (and src (list src))))))

(defun genetics--load-annotation-org (file)
  "Return the annotations in Org FILE (headings with property drawers)."
  (with-temp-buffer
    (let ((coding-system-for-read 'utf-8))
      (insert-file-contents file))
    (goto-char (point-min))
    (let ((result nil) (props nil) (body nil) (in-drawer nil) (heading nil))
      (cl-flet ((flush ()
                  (when heading
                    (let ((p props))
                      (unless (assoc "RSID" p)
                        (when (string-match "\\<\\(rs[0-9]+\\)\\>" heading)
                          (push (cons "RSID" (match-string 1 heading)) p)))
                      (when (or p (assoc "RSID" p))
                        (push (genetics--ann-from-org
                               p (string-join (nreverse body) "\n") file)
                              result))))
                  (setq props nil body nil in-drawer nil)))
        (while (not (eobp))
          (let ((line (buffer-substring-no-properties
                       (line-beginning-position) (line-end-position))))
            (cond
             ((string-match "\\`\\*+ +\\(.*\\)" line)
              (let ((title (match-string 1 line)))
                (flush)
                (setq heading title)))
             ((string-match-p "\\`[ \t]*:PROPERTIES:[ \t]*\\'" line)
              (setq in-drawer t))
             ((string-match-p "\\`[ \t]*:END:[ \t]*\\'" line)
              (setq in-drawer nil))
             ((and in-drawer
                   (string-match "\\`[ \t]*:\\([A-Za-z0-9_]+\\):[ \t]*\\(.*\\)\\'"
                                 line))
              (push (cons (upcase (match-string 1 line))
                          (string-trim (match-string 2 line)))
                    props))
             (heading (push line body))))
          (forward-line 1))
        (flush))
      (nreverse result))))

(defun genetics-load-annotation-file (file)
  "Return the list of annotations read from FILE (.json or .org)."
  (unless (file-readable-p file)
    (genetics--ann-error file "file not found or unreadable"))
  (pcase (downcase (or (file-name-extension file) ""))
    ("json" (genetics--load-annotation-json file))
    ("org" (genetics--load-annotation-org file))
    (ext (genetics--ann-error file "unsupported extension %S (use .json or .org)"
                              ext))))

(defvar genetics--annotation-cache nil
  "Cons of (SIGNATURE . HASH-TABLE) for the loaded annotations.")

(defvar genetics--annotation-position-cache nil
  "Cons of (HASH-TABLE . POSITION-TABLE) indexing annotations by site.
POSITION-TABLE maps \"BUILD:CHROM:POS\" to an annotation.")

(defun genetics--annotation-signature ()
  "Return a value that differs once annotation files are modified."
  (mapcar (lambda (f)
            (list f (file-attribute-modification-time (file-attributes f))
                  (file-attribute-size (file-attributes f))))
          genetics-annotation-files))

(defun genetics-annotations ()
  "Return a hash table of rsid to annotation for `genetics-annotation-files'."
  (let ((sig (genetics--annotation-signature)))
    (unless (and genetics--annotation-cache
                 (equal sig (car genetics--annotation-cache)))
      (let ((table (make-hash-table :test 'equal)))
        (dolist (f genetics-annotation-files)
          (dolist (a (genetics-load-annotation-file f))
            (puthash (genetics-annotation-rsid a) a table)))
        (setq genetics--annotation-cache (cons sig table))))
    (cdr genetics--annotation-cache)))

(defun genetics-reload-annotations ()
  "Discard the annotation cache so files are re-read on next use."
  (interactive)
  (setq genetics--annotation-cache nil)
  (hash-table-count (genetics-annotations)))

(defun genetics-annotation (rsid)
  "Return the annotation for RSID, or nil."
  (gethash rsid (genetics-annotations)))

(defun genetics-annotation-site (ann build)
  "Return the plist (:chrom :pos :ref) of ANN on BUILD, or nil."
  (cdr (assoc (genetics-build-key build) (genetics-annotation-coordinates ann))))

(defun genetics-annotation-at (build chrom pos)
  "Return the annotation whose BUILD coordinate is CHROM:POS, or nil."
  (let ((table (genetics-annotations)))
    (unless (eq (car genetics--annotation-position-cache) table)
      (let ((pt (make-hash-table :test 'equal)))
        (maphash (lambda (_r ann)
                   (dolist (c (genetics-annotation-coordinates ann))
                     (puthash (format "%s:%s:%d" (car c)
                                      (plist-get (cdr c) :chrom)
                                      (plist-get (cdr c) :pos))
                              ann pt)))
                 table)
        (setq genetics--annotation-position-cache (cons table pt))))
    (let ((key (genetics-build-key build)))
      (when key
        (gethash (format "%s:%s:%d" key chrom pos)
                 (cdr genetics--annotation-position-cache))))))

(defun genetics-snp-annotation (snp &optional build annotations)
  "Return the annotation of SNP by rsid, or by its position on BUILD.
ANNOTATIONS is the rsid table (default `genetics-annotations')."
  (or (gethash (genetics-snp-rsid snp) (or annotations (genetics-annotations)))
      (and build (genetics-annotation-at build (genetics-snp-chrom snp)
                                         (genetics-snp-pos snp)))))

(provide 'genetics-catalog)
;;; genetics-catalog.el ends here
