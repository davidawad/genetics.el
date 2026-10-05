;;; genetics-base.el --- Customization, errors and JSON helpers for genetics.el -*- lexical-binding: t; -*-

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

;; Options, error types, JSON and output-file helpers that every
;; genetics.el module builds on.  Nothing here touches the network.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)

;; Defined only in builds with native JSON (libjansson, on Emacs 29); each
;; use is guarded, so builds without it still compile clean.
(declare-function json-available-p "json.c" ())
(declare-function json-parse-string "json.c" (string &rest args))
(declare-function json-serialize "json.c" (object &rest args))

(defgroup genetics nil
  "Read and explore consumer genetics raw-data files."
  :group 'tools
  :prefix "genetics-")

(defconst genetics--package-dir
  (file-name-directory (or load-file-name buffer-file-name default-directory))
  "Directory the genetics package was loaded from.")

;;;; Customization

(defconst genetics--conventional-data-directory "~/Documents/Genetics/"
  "Conventional folder for raw-data files, the default when it exists.")

(defun genetics--default-data-directory ()
  "Return the default for `genetics-data-directory' on this system.
That is `genetics--conventional-data-directory' when it exists, else
the home directory.  Set `genetics-data-directory' to point elsewhere."
  (if (file-directory-p genetics--conventional-data-directory)
      genetics--conventional-data-directory
    "~/"))

(defcustom genetics-data-directory (genetics--default-data-directory)
  "Default directory offered when prompting for a raw-data file.
When it does not exist (a configuration shared between machines, say)
prompts start in `default-directory' instead."
  :type 'directory)

(defun genetics--prompt-directory ()
  "Return `genetics-data-directory' if it exists, else `default-directory'."
  (if (and (stringp genetics-data-directory)
           (file-directory-p genetics-data-directory))
      (file-name-as-directory genetics-data-directory)
    default-directory))

(defcustom genetics-cache-directory (locate-user-emacs-file "genetics-cache/")
  "Directory for parse caches, decompressed VCFs and SNPedia answers."
  :type 'directory)

(defcustom genetics-gzip-program "gzip"
  "Name or path of the gzip executable used to read .gz files.
It is looked up with `executable-find'.  When it is nil or not found,
.gz files are decompressed with Emacs' built-in zlib support instead
\(see `zlib-available-p'), which reads the whole file into memory."
  :type '(choice (string :tag "Executable")
                 (const :tag "Always use Emacs' zlib" nil)))

(defcustom genetics-use-cache t
  "Non-nil means cache parsed kits on disk as `.eld' files."
  :type 'boolean)

(defcustom genetics-vcf-eager-limit (* 200 1024 1024)
  "VCF files larger than this many bytes are offset-indexed, not loaded."
  :type 'integer)

(defcustom genetics-chunk-size (* 4 1024 1024)
  "Bytes read at a time when streaming a large file."
  :type 'integer)

(defcustom genetics-browse-limit 5000
  "Maximum number of rows shown in the browser buffer."
  :type 'integer)

(defcustom genetics-annotation-files
  (list (expand-file-name "annotations/genetics-curated.json"
                          genetics--package-dir))
  "List of annotation files (JSON or Org); later files override earlier."
  :type '(repeat file))

(defcustom genetics-vcf-ref-calls 'auto
  "How absent sites of a VCF parsed in Emacs are treated.
`auto' treats a VCF as variant-only whole-genome data (`absent-means-ref')
when it has no explicit homozygous-reference calls and at least
`genetics-wgs-min-records' records, and as `explicit' when it is a gVCF;
otherwise `unknown'.  `absent-means-ref' and `unknown' force that
answer.  Only `absent-means-ref' lets curated sites that are missing from
the file be shown as inferred homozygous reference."
  :type '(choice (const :tag "Detect" auto)
                 (const :tag "Variant-only WGS (absent = reference)"
                        absent-means-ref)
                 (const :tag "Never infer" unknown)))

(defcustom genetics-wgs-min-records 1000000
  "Minimum VCF record count for `auto' to treat a VCF as whole-genome."
  :type 'integer)

(defcustom genetics-snpedia-enabled nil
  "Non-nil allows SNPedia lookups, which send the rsid over the network."
  :type 'boolean)

;;;; Errors

(define-error 'genetics-error "Genetics error")
(define-error 'genetics-file-error "Genetics file error" 'genetics-error)
(define-error 'genetics-unknown-format "Unrecognized genetics file format"
  'genetics-error)
(define-error 'genetics-parse-error "Cannot parse genetics data"
  'genetics-error)
(define-error 'genetics-gzip-error "Cannot decompress gzip file"
  'genetics-error)
(define-error 'genetics-unsupported-file
  "Raw reads cannot be opened as genotypes" 'genetics-error)
(define-error 'genetics-fastq-file "FASTQ reads cannot be opened directly"
  'genetics-unsupported-file)
(define-error 'genetics-genome-error "genome-cli failed" 'genetics-error)
(define-error 'genetics-genome-missing "genome-cli executable not found"
  'genetics-genome-error)
(define-error 'genetics-no-kit "No genetics kit available" 'genetics-error)
(define-error 'genetics-annotation-error "Bad genetics annotation"
  'genetics-error)
(define-error 'genetics-snpedia-error "SNPedia lookup failed" 'genetics-error)
(define-error 'genetics-snpedia-disabled
  "SNPedia lookups are disabled (set `genetics-snpedia-enabled')"
  'genetics-snpedia-error)
(define-error 'genetics-snpedia-declined "SNPedia lookup not confirmed"
  'genetics-snpedia-error)

(defun genetics--native-json-p ()
  "Return non-nil if Emacs has native JSON (always on 30, libjansson on 29)."
  (and (fboundp 'json-available-p) (json-available-p)))

(defun genetics--json-parse (string &rest args)
  "Parse JSON STRING like `json-parse-string' with keyword ARGS.
Uses native JSON when available, else json.el, so Emacs 29 builds
without libjansson work too.  Errors are `json-error' either way."
  (if (genetics--native-json-p)
      (apply #'json-parse-string string args)
    (let ((json-object-type (pcase (plist-get args :object-type)
                              ('alist 'alist) ('plist 'plist) (_ 'hash-table)))
          (json-array-type (if (eq (plist-get args :array-type) 'list) 'list 'vector))
          (json-key-type nil)
          (json-null (if (plist-member args :null-object)
                         (plist-get args :null-object)
                       :null))
          (json-false (if (plist-member args :false-object)
                          (plist-get args :false-object)
                        :false)))
      (json-read-from-string string))))

(defun genetics--json-serialize (object)
  "Return OBJECT as a JSON string like `json-serialize' (:null, :false)."
  (if (genetics--native-json-p)
      (json-serialize object)
    (let ((json-null :null) (json-false :false)
          (json-encoding-pretty-print nil))
      (json-encode object))))

(defmacro genetics--with-output-file (file &rest body)
  "Like `with-temp-file' on FILE with BODY, but always write UTF-8 with LF.
Data files are byte-for-byte the same on every OS (no CRLF, no locale
coding on Windows)."
  (declare (indent 1) (debug t))
  `(with-temp-file ,file
     (setq buffer-file-coding-system 'utf-8-unix)
     ,@body))

(defun genetics--error (type fmt &rest args)
  "Signal error TYPE with a message built from FMT and ARGS."
  (signal type (list (apply #'format fmt args))))

(provide 'genetics-base)
;;; genetics-base.el ends here
