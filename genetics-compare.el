;;; genetics-compare.el --- Compare two genetics kits -*- lexical-binding: t; -*-

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

;; `genetics-compare' reports overlap and concordance between two kits,
;; ignoring no-calls and allele order, and flags complement-strand
;; discordances separately.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'genetics-core)
(require 'genetics-stats)

(defun genetics--complement-key (key)
  "Return the genotype key KEY with every allele complemented."
  (genetics-genotype-key
   (mapconcat (lambda (a) (or (genetics-complement a) a))
              (genetics-alleles key) "/")))

(defun genetics-compare-kits (a b)
  "Compare kits A and B; return a plist of results.
Keys: :overlap :compared :concordant :discordant (list of (RSID GA GB))
:complement (list of (RSID GA GB)) :concordance (percent or nil)
:build-mismatch.  Only records with calls in both kits are compared;
hemizygous alleles compare as their homozygous form.  Two kits served
by the same genome-cli backend are compared by genome-cli."
  (let ((a (genetics-find-kit a)) (b (genetics-find-kit b)))
    (cond
     ((and (genetics-kit-backend a)
           (eq (genetics-kit-backend a) (genetics-kit-backend b)))
      (funcall (genetics-kit-backend a) 'compare a b))
     ((or (genetics-kit-backend a) (genetics-kit-backend b))
      (genetics--error 'genetics-error
                       "Cannot compare a genome-cli kit with a kit parsed in Emacs; open both with the same `genetics-source-function'"))
     (t (genetics--compare-native a b)))))

(defun genetics--compare-native (a b)
  "Compare kits A and B parsed in Emacs, see `genetics-compare-kits'."
  (let* ((small (if (<= (genetics-kit-snp-count a) (genetics-kit-snp-count b))
                    a b))
         (overlap 0) (compared 0) (concordant 0) (discordant nil)
         (complement nil))
    (genetics-kit-map-snps
     small
     (lambda (snp)
       (let* ((id (genetics-snp-rsid snp))
              (sa (if (eq small a) snp (genetics-kit-get a id)))
              (sb (if (eq small b) snp (genetics-kit-get b id))))
         (when (and sa sb)
           (cl-incf overlap)
           (let ((ka (genetics-genotype-key (genetics-snp-genotype sa)))
                 (kb (genetics-genotype-key (genetics-snp-genotype sb))))
             (when (and ka kb)
               (cl-incf compared)
               (cond ((equal ka kb) (cl-incf concordant))
                     ((equal ka (genetics--complement-key kb))
                      (push (list id (genetics-snp-genotype sa)
                                  (genetics-snp-genotype sb))
                            complement))
                     (t (push (list id (genetics-snp-genotype sa)
                                    (genetics-snp-genotype sb))
                              discordant)))))))
       nil))
    (list :overlap overlap :compared compared :concordant concordant
          :discordant (nreverse discordant) :complement (nreverse complement)
          :concordance (when (> compared 0)
                         (* 100.0 (/ (float concordant) compared)))
          :build-mismatch (genetics-builds-differ-p a b))))

;;;###autoload
(defun genetics-compare (kit-a kit-b)
  "Compare KIT-A and KIT-B and show the result in a buffer.
Returns the buffer."
  (interactive
   (progn
     (when (< (length genetics-loaded-kits) 2)
       (genetics--error 'genetics-no-kit "Need two loaded kits to compare"))
     (let ((names (mapcar #'genetics-kit-name genetics-loaded-kits)))
       (list (completing-read "Kit A: " names nil t)
             (completing-read "Kit B: " names nil t)))))
  (let* ((a (genetics-find-kit kit-a)) (b (genetics-find-kit kit-b))
         (r (genetics-compare-kits a b))
         (buf (get-buffer-create (format "*genetics-compare: %s vs %s*"
                                         (genetics-kit-name a)
                                         (genetics-kit-name b)))))
    (with-current-buffer buf
      (special-mode)
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert (format "Comparing %s (GRCh%s) with %s (GRCh%s)\n\n"
                        (genetics-kit-name a) (or (genetics-kit-build a) "?")
                        (genetics-kit-name b) (or (genetics-kit-build b) "?")))
        (when (plist-get r :build)
          (insert (format "Compared by genome-cli on %s.\n" (plist-get r :build))))
        (dolist (w (plist-get r :warnings))
          (insert (format "Note: %s\n" w)))
        (when (or (plist-get r :build) (plist-get r :warnings))
          (insert "\n"))
        (when (plist-get r :build-mismatch)
          (insert "WARNING: the kits are on different genome builds. Matching is by rsid or chrom:pos id; positions and reference alleles differ between builds and no liftover is performed, so results may be misleading.\n\n"))
        (insert (format "Overlapping ids:     %d\n" (plist-get r :overlap)))
        (insert (format "Compared (both called): %d\n" (plist-get r :compared)))
        (insert (format "Concordant:          %d\n" (plist-get r :concordant)))
        (insert (format "Concordance:         %s\n"
                        (if (plist-get r :concordance)
                            (format "%.2f%%" (plist-get r :concordance))
                          "n/a")))
        (insert (format "Complement-strand:   %d (same alleles on the opposite strand; excluded from concordant)\n"
                        (length (plist-get r :complement))))
        (insert (format "Discordant:          %d\n\n"
                        (or (plist-get r :discordant-count)
                            (length (plist-get r :discordant)))))
        (dolist (section '((:discordant . "Discordant calls")
                           (:complement . "Possible complement-strand differences")))
          (when (plist-get r (car section))
            (insert (format "%s\n" (cdr section)))
            (dolist (d (plist-get r (car section)))
              (insert (format "  %-14s %-6s %-6s\n" (nth 0 d) (nth 1 d) (nth 2 d))))
            (insert "\n")))
        (goto-char (point-min))))
    (pop-to-buffer buf)
    buf))

(provide 'genetics-compare)
;;; genetics-compare.el ends here
