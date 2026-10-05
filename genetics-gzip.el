;;; genetics-gzip.el --- Gzip and BGZF decompression for genetics.el -*- lexical-binding: t; -*-

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

;; Decompresses .gz and bgzip (BGZF) files with the local gzip executable
;; or Emacs' own zlib, and keys cached copies by file identity.  Nothing
;; is sent anywhere.

;;; Code:

(require 'genetics-core)

;;;; gzip

;; Defined only in builds with zlib; guarded by `genetics--zlib-p'.
(declare-function zlib-available-p "decompress.c" ())
(declare-function zlib-decompress-region "decompress.c"
                  (start end &optional allow-partial))

(defun genetics--gzip ()
  "Return the gzip executable, or nil if it is not installed."
  (and genetics-gzip-program (executable-find genetics-gzip-program)))

(defun genetics--zlib-p ()
  "Return non-nil if this Emacs can decompress gzip data itself."
  (and (fboundp 'zlib-available-p) (zlib-available-p)))

(defun genetics--no-gunzip-error (file)
  "Signal that FILE cannot be decompressed here."
  (genetics--error 'genetics-gzip-error
                   "Cannot read %s: no gzip executable and this Emacs lacks zlib; install gzip, use genome-cli, or decompress the file first"
                   file))

(defun genetics--u16 (pos)
  "Return the little-endian 16-bit integer at POS in the current buffer."
  (+ (char-after pos) (* 256 (char-after (1+ pos)))))

(defun genetics--bgzf-size (start)
  "Return the size of the BGZF block at START, or nil if not BGZF.
BGZF (bgzip, tabix) members carry their length in a \"BC\" extra field."
  (when (and (<= (+ start 12) (point-max))
             (/= 0 (logand 4 (char-after (+ start 3)))))
    (let* ((p (+ start 12))
           (xend (min (point-max) (+ p (genetics--u16 (+ start 10)))))
           (size nil))
      (while (and (not size) (<= (+ p 4) xend))
        (let ((slen (genetics--u16 (+ p 2))))
          (if (and (= (char-after p) ?B) (= (char-after (1+ p)) ?C) (= slen 2)
                   (<= (+ p 6) xend))
              (setq size (1+ (genetics--u16 (+ p 4))))
            (setq p (+ p 4 slen)))))
      size)))

(defun genetics--zlib-gunzip-into-buffer (file)
  "Insert the decompressed contents of gzip FILE at point, using zlib.
Handles single-member gzip and multi-member BGZF; signals
`genetics-gzip-error' for other multi-member files or corrupt data."
  (let ((out (current-buffer)) (start 1))
    (with-temp-buffer
      (set-buffer-multibyte nil)
      (insert-file-contents-literally file)
      (let ((src (current-buffer)))
        (while (< start (point-max))
          (unless (and (<= (+ start 10) (point-max))
                       (= (char-after start) #x1f)
                       (= (char-after (1+ start)) #x8b))
            (genetics--error 'genetics-gzip-error "Not gzip data in %s" file))
          (let* ((size (genetics--bgzf-size start))
                 (end (if size (min (point-max) (+ start size)) (point-max)))
                 (isize (and (null size) (>= (- end start) 18)
                             (+ (genetics--u16 (- end 4))
                                (* 65536 (genetics--u16 (- end 2)))))))
            (with-temp-buffer
              (set-buffer-multibyte nil)
              (insert-buffer-substring src start end)
              (unless (zlib-decompress-region (point-min) (point-max))
                (genetics--error 'genetics-gzip-error
                                 "Corrupt gzip data in %s" file))
              (when (and isize (/= isize (logand (buffer-size) #xFFFFFFFF)))
                (genetics--error 'genetics-gzip-error
                                 "%s has several gzip members; install gzip to read it"
                                 file))
              (let ((chunk (current-buffer)))
                (with-current-buffer out (insert-buffer-substring chunk))))
            (setq start end)))))))

(defun genetics--gunzip-into-buffer (file)
  "Insert the decompressed contents of FILE at point.
Uses `genetics-gzip-program' when installed, else Emacs' zlib."
  (let ((gzip (genetics--gzip)))
    (cond
     (gzip
      (let* ((coding-system-for-read 'no-conversion)
             (status (call-process gzip nil t nil "-dc" "--" file)))
        (unless (eq status 0)
          (genetics--error 'genetics-gzip-error
                           "gzip failed on %s (exit status %s)" file status))))
     ((genetics--zlib-p) (genetics--zlib-gunzip-into-buffer file))
     (t (genetics--no-gunzip-error file)))))

(defun genetics--cache-key (file)
  "Return a hash string identifying FILE by truename, size and mtime."
  (let* ((truename (file-truename file))
         (attrs (file-attributes truename)))
    (secure-hash 'sha1 (format "%s|%d|%s" truename
                               (file-attribute-size attrs)
                               (format-time-string
                                "%s.%N" (file-attribute-modification-time
                                         attrs))))))

(defun genetics--decompressed-copy (file)
  "Return the path of a decompressed copy of gz FILE in the cache directory."
  (let ((out (expand-file-name (format "vcf-%s.vcf" (genetics--cache-key file))
                               genetics-cache-directory)))
    (unless (and (file-exists-p out)
                 (> (file-attribute-size (file-attributes out)) 0))
      (make-directory genetics-cache-directory t)
      (let ((tmp (concat out ".part"))
            (gzip (genetics--gzip)))
        (cond
         (gzip
          (let ((status (call-process gzip nil (list :file tmp) nil
                                      "-dc" "--" file)))
            (unless (eq status 0)
              (ignore-errors (delete-file tmp))
              (genetics--error 'genetics-gzip-error
                               "gzip failed on %s (exit status %s)" file status))))
         ((genetics--zlib-p)
          (with-temp-buffer
            (set-buffer-multibyte nil)
            (genetics--zlib-gunzip-into-buffer file)
            (let ((coding-system-for-write 'no-conversion))
              (write-region nil nil tmp nil 'silent))))
         (t (genetics--no-gunzip-error file)))
        (rename-file tmp out t)))
    out))

(provide 'genetics-gzip)
;;; genetics-gzip.el ends here
