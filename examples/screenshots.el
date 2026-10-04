;;; screenshots.el --- capture the README screenshots -*- lexical-binding: t; -*-

;;; Commentary:

;; Loaded by examples/screenshots.sh into a graphical `emacs -Q' running
;; on an X display (Xvfb).  Opens the synthetic fixtures only, shows each
;; view full-frame and captures the Emacs window with `xwd', converted to
;; PNG with ImageMagick, into docs/screenshots/.
;;
;; Environment: GENETICS_SHOT_DIR (output directory), GENETICS_SHOT_THEME
;; (default modus-vivendi-tinted), GENETICS_SHOT_FONT (default the first
;; installed of JetBrains Mono, DejaVu Sans Mono, monospace).

;;; Code:

(defvar shot-root
  (file-name-directory
   (directory-file-name (file-name-directory (or load-file-name buffer-file-name))))
  "Repository root.")

(defvar shot-dir
  (file-name-as-directory
   (or (getenv "GENETICS_SHOT_DIR") (expand-file-name "docs/screenshots" shot-root)))
  "Where PNGs are written.")

(defun shot-fixture (name)
  "Return the path of synthetic fixture NAME."
  (expand-file-name (concat "test/fixtures/" name) shot-root))

;;;; Frame look

(menu-bar-mode -1)
(tool-bar-mode -1)
(scroll-bar-mode -1)
(blink-cursor-mode -1)
(setq inhibit-startup-screen t
      ring-bell-function #'ignore
      frame-resize-pixelwise t
      use-dialog-box nil
      ;; The fixtures live in the repository; show them as a user's files.
      directory-abbrev-alist (list (cons (concat "\\`" (regexp-quote shot-root))
                                         "~/src/genetics-el/")))
(let ((font (or (getenv "GENETICS_SHOT_FONT")
                (seq-find (lambda (f) (member f (font-family-list)))
                          '("JetBrains Mono" "DejaVu Sans Mono"))
                "monospace")))
  (set-face-attribute 'default nil :family font :height 110)
  (set-face-attribute 'fixed-pitch nil :family font)
  (set-face-attribute 'variable-pitch nil :family font))
(load-theme (intern (or (getenv "GENETICS_SHOT_THEME") "modus-vivendi-tinted")) t)
(set-frame-parameter nil 'internal-border-width 14)
(setq-default mode-line-format
              '(" " mode-line-buffer-identification "   " mode-name))

(defun shot-frame-size (cols rows)
  "Resize the frame to COLS by ROWS characters and let X settle."
  (set-frame-size nil cols rows)
  (redisplay t)
  (sit-for 0.4))

;;;; Capture

(defun shot-capture (name)
  "Capture the Emacs frame as NAME.png in `shot-dir'."
  (message nil)
  (redisplay t)
  (sit-for 0.6)
  (redisplay t)
  (let* ((xwd (make-temp-file "genetics-shot" nil ".xwd"))
         (png (expand-file-name (concat name ".png") shot-dir))
         (convert (or (executable-find "magick") (executable-find "convert"))))
    (unless (zerop (call-process "xwd" nil nil nil "-silent" "-id"
                                 (frame-parameter nil 'outer-window-id)
                                 "-out" xwd))
      (error "Xwd failed for %s" name))
    (unless (zerop (call-process convert nil nil nil xwd "-strip" png))
      (error "ImageMagick failed for %s" name))
    (delete-file xwd)
    (princ (format "wrote %s\n" png) #'external-debugging-output)))

(defun shot-show (buffer)
  "Show BUFFER alone in the frame, point at the top."
  (delete-other-windows)
  (switch-to-buffer buffer)
  (setq-local cursor-type nil)
  (setq-local word-wrap t)
  (setq-local fringe-indicator-alist
              (cons '(continuation nil nil) fringe-indicator-alist))
  (goto-char (point-min)))

;;;; The views

(add-to-list 'load-path shot-root)
(require 'genetics)
(setq genetics-source-function #'genetics-source-native
      genetics-use-cache nil)
(make-directory shot-dir t)

(condition-case err
    (let* ((kit (genetics-open (shot-fixture "23andme-sample.txt"))))
      (genetics-open (shot-fixture "myheritage-sample.csv"))
      (genetics-open (shot-fixture "ancestry-sample.txt"))

      ;; 1. Summary buffer
      (shot-frame-size 96 36)
      (shot-show (genetics-summary kit))
      (shot-capture "summary")

      ;; 2. Browser, filtered to annotated SNPs
      (shot-frame-size 96 13)
      (shot-show (genetics-browse kit))
      (genetics-browse-filter-annotated)
      (goto-char (point-min))
      (shot-capture "browse")

      ;; 3. Lookup of rs429358 across the three loaded kits
      (shot-frame-size 96 36)
      (shot-show (genetics-lookup "rs429358"))
      (shot-capture "lookup")

      ;; 4. The Org report built from dynamic blocks, from its first heading
      (shot-frame-size 158 54)
      (shot-show (let ((directory-abbrev-alist nil))
                   (find-file-noselect
                    (expand-file-name "examples/genetics-report.org" shot-root))))
      (setq-local org-hide-emphasis-markers t)
      (setq-local word-wrap nil)
      (setq truncate-lines t)
      (font-lock-flush)
      (org-table-map-tables #'org-table-shrink t)
      (re-search-forward "^\\* Kit summary")
      (set-window-start nil (line-beginning-position))
      (shot-capture "org-report")
      (kill-emacs 0))
  (error (princ (format "screenshot failed: %S\n" err) #'external-debugging-output)
         (kill-emacs 1)))

;;; screenshots.el ends here
