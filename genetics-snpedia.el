;;; genetics-snpedia.el --- Opt-in SNPedia lookups for genetics.el -*- lexical-binding: t; -*-

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

;; PRIVACY: this is the ONLY file in genetics.el that talks to the network.
;; Lookups are disabled unless `genetics-snpedia-enabled' is non-nil, ask for
;; confirmation on first use, and send ONLY the rsid (a public identifier,
;; never genotypes, positions or file names).  Answers are cached on disk.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'json)
(require 'url)
(require 'genetics-core)

(defcustom genetics-snpedia-url-format
  "https://bots.snpedia.com/api.php?action=parse&page=%s&format=json&prop=wikitext"
  "Request URL; the single %s is replaced by the capitalized rsid."
  :type 'string
  :group 'genetics)

(defcustom genetics-snpedia-timeout 15
  "Seconds to wait for a SNPedia response."
  :type 'integer
  :group 'genetics)

(defvar genetics-snpedia-confirmed nil
  "Non-nil once the user confirmed SNPedia lookups in this session.")

(defun genetics-snpedia-url (rsid)
  "Return the request URL for RSID, which must look like rs123."
  (unless (and (stringp rsid) (string-match-p "\\`rs[0-9]+\\'" rsid))
    (genetics--error 'genetics-snpedia-error
                     "Only rsids (rs followed by digits) can be looked up, got %S"
                     rsid))
  (format genetics-snpedia-url-format (concat "R" (substring rsid 1))))

(defun genetics--snpedia-cache-file (rsid)
  "Return the cache file for RSID."
  (expand-file-name (concat rsid ".json")
                    (expand-file-name "snpedia/" genetics-cache-directory)))

(defun genetics--snpedia-confirm ()
  "Ask once per session for permission to contact SNPedia."
  (unless genetics-snpedia-confirmed
    (if (y-or-n-p "Look up this rsid on SNPedia (sends only the rsid to bots.snpedia.com)? ")
        (setq genetics-snpedia-confirmed t)
      (genetics--error 'genetics-snpedia-declined
                       "SNPedia lookup cancelled; nothing was sent"))))

(defun genetics--snpedia-download (rsid)
  "Fetch the raw JSON text for RSID from SNPedia."
  (let* ((url (genetics-snpedia-url rsid))
         (url-privacy-level 'paranoid)
         (url-automatic-caching nil)
         (buf (url-retrieve-synchronously url t t genetics-snpedia-timeout)))
    (unless buf
      (genetics--error 'genetics-snpedia-error
                       "No response from SNPedia for %s" rsid))
    (unwind-protect
        (with-current-buffer buf
          (goto-char (point-min))
          (unless (looking-at "HTTP/[0-9.]+ 200")
            (genetics--error 'genetics-snpedia-error
                             "SNPedia returned an error for %s: %s" rsid
                             (buffer-substring (point-min) (line-end-position))))
          (unless (re-search-forward "\r?\n\r?\n" nil t)
            (genetics--error 'genetics-snpedia-error
                             "Malformed SNPedia response for %s" rsid))
          (decode-coding-string
           (buffer-substring-no-properties (point) (point-max)) 'utf-8))
      (kill-buffer buf))))

(defun genetics--snpedia-wikitext (json-text rsid)
  "Extract the wikitext string for RSID from JSON-TEXT."
  (let* ((data (condition-case nil
                   (genetics--json-parse json-text :object-type 'alist)
                 (error (genetics--error 'genetics-snpedia-error
                                         "SNPedia sent invalid JSON for %s" rsid))))
         (parse (alist-get 'parse data))
         (wt (alist-get 'wikitext parse)))
    (cond ((stringp wt) wt)
          ((stringp (alist-get '* wt)) (alist-get '* wt))
          (t (genetics--error 'genetics-snpedia-error
                              "SNPedia has no page for %s" rsid)))))

(defun genetics-snpedia-fetch (rsid)
  "Return the SNPedia wikitext for RSID (cached on disk).
Signals `genetics-snpedia-disabled' unless `genetics-snpedia-enabled'."
  (unless genetics-snpedia-enabled
    (signal 'genetics-snpedia-disabled nil))
  (let ((cache (genetics--snpedia-cache-file rsid)))
    (genetics-snpedia-url rsid)         ; validate before anything else
    (if (file-readable-p cache)
        (genetics--snpedia-wikitext
         (with-temp-buffer
           (let ((coding-system-for-read 'utf-8))
             (insert-file-contents cache))
           (buffer-string)) rsid)
      (genetics--snpedia-confirm)
      (let ((text (genetics--snpedia-download rsid)))
        (let ((wikitext (genetics--snpedia-wikitext text rsid)))
          (make-directory (file-name-directory cache) t)
          (genetics--with-output-file cache (insert text))
          wikitext)))))

(defun genetics--snpedia-field (wikitext name)
  "Return template field NAME from WIKITEXT, or nil."
  (let ((case-fold-search t))
    (when (string-match (format "|[ \t]*%s[ \t]*=[ \t]*\\([^|\n}]*\\)" name)
                        wikitext)
      (let ((v (string-trim (match-string 1 wikitext))))
        (unless (string-empty-p v) v)))))

;;;###autoload
(defun genetics-snpedia-summary (rsid)
  "Return a plist (:rsid :gene :orientation :summary :wikitext) for RSID."
  (let ((wt (genetics-snpedia-fetch rsid)))
    (list :rsid rsid
          :gene (genetics--snpedia-field wt "gene")
          :orientation (genetics--snpedia-field wt "orientation")
          :summary (genetics--snpedia-field wt "summary")
          :wikitext wt)))

(provide 'genetics-snpedia)
;;; genetics-snpedia.el ends here
