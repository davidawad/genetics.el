;;; regenerate-report.el --- rebuild examples/genetics-report.{org,html} -*- lexical-binding: t; -*-

;;; Commentary:

;; Run from the repository root:
;;
;;   emacs -Q --batch -l examples/regenerate-report.el
;;
;; Updates every dynamic block of examples/genetics-report.org from the
;; synthetic fixture test/fixtures/23andme-sample.txt with the Emacs Lisp
;; parser (so the output does not depend on genome-cli being installed),
;; saves it, and exports it to examples/genetics-report.html.

;;; Code:

(let* ((root (file-name-directory
              (directory-file-name
               (file-name-directory (or load-file-name buffer-file-name)))))
       (org-file (expand-file-name "examples/genetics-report.org" root)))
  (add-to-list 'load-path root)
  (require 'genetics)
  (require 'ox-html)
  (setq genetics-source-function #'genetics-source-native
        genetics-use-cache nil
        genetics-loaded-kits nil)
  (with-current-buffer (find-file-noselect org-file)
    ;; Keep the machine's own paths out of the sample.
    (let ((directory-abbrev-alist
           (list (cons (concat "\\`" (regexp-quote root)) "~/src/genetics-el/"))))
      (org-update-all-dblocks))
    (save-buffer)
    (let ((org-html-validation-link nil)
          (org-export-time-stamp-file nil)
          (org-html-postamble nil))
      ;; Org numbers headings with `random'; seed it for stable diffs.
      (random "genetics-report")
      (org-html-export-to-html))
    (message "Wrote %s and its HTML export" (abbreviate-file-name org-file))))

;;; regenerate-report.el ends here
