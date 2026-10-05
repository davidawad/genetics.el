;;; genetics-fastq.el --- FASTQ reads to a VCF through genome-cli -*- lexical-binding: t; -*-

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

;; genetics.el does not parse FASTQ: reads must be aligned and
;; variant-called first.  These commands drive genome-cli's pipeline:
;;
;;   M-x genetics-fastq-plan   `genome pipeline plan', shown as a table
;;   M-x genetics-fastq-run    show the plan, confirm, then run
;;                             `genome pipeline run' asynchronously in a
;;                             compilation buffer; on success offer to
;;                             open the resulting VCF
;;
;; `genetics-fastq-plan-explain' and `genetics-fastq-run-explain' return
;; the exact command lines without running anything.  The pipeline runs
;; locally; its fetch-reference step (if not cached) downloads the public
;; reference genome, which sends none of your data.

;;; Code:

(require 'cl-lib)
(require 'subr-x)
(require 'compile)
(require 'json)
(require 'genetics-core)
(require 'genetics-source)

(declare-function genetics-open "genetics" (file))

(defcustom genetics-fastq-build "GRCh38"
  "Reference build the FASTQ pipeline aligns to."
  :type '(choice (const "GRCh38") (const "GRCh37"))
  :group 'genetics)

(defcustom genetics-fastq-output-directory nil
  "Directory for pipeline outputs, or nil for genome-cli's default."
  :type '(choice (const :tag "genome-cli default" nil) directory)
  :group 'genetics)

(defcustom genetics-fastq-extra-args nil
  "Extra arguments passed to `genome pipeline plan' and `run'."
  :type '(repeat string)
  :group 'genetics)

;;;; argv (pure)

(defun genetics-fastq-argv (action files)
  "Return the genome-cli argv for pipeline ACTION (`plan' or `run') on FILES.
FILES are FASTQ paths (one file, or R1 and R2 of a pair).  Pure."
  (unless (memq action '(plan run))
    (genetics--error 'genetics-genome-error "Unknown pipeline action %S" action))
  (unless files
    (genetics--error 'genetics-file-error "No FASTQ files given"))
  (append (list genetics-genome-executable "pipeline" (symbol-name action))
          (mapcar #'expand-file-name files)
          (list "--build" genetics-fastq-build)
          (when genetics-fastq-output-directory
            (list "--out" (expand-file-name genetics-fastq-output-directory)))
          genetics-fastq-extra-args
          (list "--format" "json")))

(defun genetics--read-fastq-files ()
  "Read one or two FASTQ file names (R1 and optional R2)."
  (let* ((r1 (read-file-name "FASTQ (R1): " (genetics--prompt-directory) nil t))
         (r2 (read-file-name "Mate FASTQ (R2, empty if single-end): "
                             (file-name-directory r1) "" nil)))
    (if (or (string-empty-p r2) (genetics--same-file-p r2 r1)
            (file-directory-p r2))
        (list r1)
      (list r1 r2))))

(defun genetics-fastq-plan-explain (files)
  "Return the exact `genome pipeline plan' command for FILES.
Pure twin of `genetics-fastq-plan'; interactively, show it."
  (interactive (list (genetics--read-fastq-files)))
  (let ((cmd (genetics-genome-argv-string (genetics-fastq-argv 'plan files))))
    (when (called-interactively-p 'interactive) (message "%s" cmd))
    cmd))

(defun genetics-fastq-run-explain (files)
  "Return the exact `genome pipeline run' command for FILES.
Pure twin of `genetics-fastq-run'; interactively, show it."
  (interactive (list (genetics--read-fastq-files)))
  (let ((cmd (genetics-genome-argv-string (genetics-fastq-argv 'run files))))
    (when (called-interactively-p 'interactive) (message "%s" cmd))
    cmd))

;;;; Plan buffer

(defvar-local genetics-fastq--files nil
  "FASTQ files of the current plan buffer.")

(defvar-local genetics-fastq--steps nil
  "Pipeline steps (genome/v1 alists) of the current plan buffer.")

(defvar genetics-fastq-plan-mode-map
  (let ((map (make-sparse-keymap)))
    (set-keymap-parent map special-mode-map)
    (define-key map "x" #'genetics-fastq-plan-run)
    map)
  "Keymap for `genetics-fastq-plan-mode'.")

(define-derived-mode genetics-fastq-plan-mode special-mode "Genetics-FASTQ-Plan"
  "Major mode showing a genome-cli FASTQ pipeline plan.

\\{genetics-fastq-plan-mode-map}")

(defun genetics-fastq--vcf-output (step-list)
  "Return the last VCF output named in STEP-LIST, or nil."
  (let (vcf)
    (dolist (s step-list)
      (dolist (o (alist-get 'outputs s))
        (when (and (stringp o) (string-match-p "\\.g?vcf\\(\\.gz\\)?\\'" o))
          (setq vcf o))))
    vcf))

(defun genetics-fastq--insert-steps (step-list)
  "Insert a description of each pipeline step in STEP-LIST."
  (let ((i 0))
    (dolist (s step-list)
      (cl-incf i)
      (insert (format "%d. %-15s %-12s %s%s\n" i
                      (alist-get 'step s) (or (alist-get 'tool s) "")
                      (or (alist-get 'status s) "")
                      (let ((sec (alist-get 'seconds s)))
                        (if (numberp sec) (format " (%.0f s)" sec) ""))))
      (when (alist-get 'argv s)
        (insert "   $ " (genetics-genome-argv-string (alist-get 'argv s)) "\n"))
      (dolist (k '((inputs . "in ") (outputs . "out")))
        (dolist (f (alist-get (car k) s))
          (insert (format "   %s %s\n" (cdr k) f)))))))

(defun genetics-fastq-plan (files)
  "Show genome-cli's plan for turning FASTQ FILES into a VCF.
FILES is a list of one FASTQ file or an R1/R2 pair.  Nothing is run
except `genome pipeline plan'.  In the plan buffer,
\\<genetics-fastq-plan-mode-map>\\[genetics-fastq-plan-run] runs it.  Returns the buffer."
  (interactive (list (genetics--read-fastq-files)))
  (let* ((argv (genetics-fastq-argv 'plan files))
         (env (genetics-genome-run argv "pipeline-plan"))
         (steps (alist-get 'data env))
         (buf (get-buffer-create "*genetics-fastq-plan*")))
    (with-current-buffer buf
      (genetics-fastq-plan-mode)
      (setq genetics-fastq--files (mapcar #'expand-file-name files)
            genetics-fastq--steps steps)
      (let ((inhibit-read-only t))
        (erase-buffer)
        (insert (propertize "FASTQ -> VCF pipeline plan (genome-cli)\n\n" 'face 'bold))
        (insert "Reads:     " (string-join (mapcar #'abbreviate-file-name files) "\n           ") "\n")
        (insert "Build:     " genetics-fastq-build "\n")
        (insert "Planned:   $ " (genetics-genome-argv-string argv) "\n")
        (insert "Run with:  $ " (genetics-fastq-run-explain files) "\n\n")
        (genetics-fastq--insert-steps steps)
        (let ((vcf (genetics-fastq--vcf-output steps)))
          (insert "\nResulting VCF: " (or vcf "(not stated in the plan)") "\n"))
        (dolist (w (alist-get 'warnings env))
          (insert "Warning: " w "\n"))
        (insert "\nEverything runs on this machine. A fetch-reference step that is not cached downloads the public reference genome; none of your data is sent. Alignment and variant calling of a 30x genome take hours and need tens of GB of disk.\n\n")
        (insert (substitute-command-keys
                 "Press \\<genetics-fastq-plan-mode-map>\\[genetics-fastq-plan-run] to run this pipeline.\n"))
        (goto-char (point-min))))
    (pop-to-buffer buf)
    buf))

;;;; Running

(defvar-local genetics-fastq--planned-vcf nil
  "VCF the plan of the current run buffer said it would produce.")

(define-compilation-mode genetics-fastq-run-mode "Genetics-FASTQ-Run"
  "Compilation mode for a running `genome pipeline run'.")

(defun genetics-fastq--result-steps (buffer)
  "Return the data of the last pipeline-run envelope printed in BUFFER."
  (with-current-buffer buffer
    (save-excursion
      (goto-char (point-max))
      (let (steps)
        (while (and (not steps)
                    (re-search-backward "^{\"schema\":\"genome/v1\"" nil t))
          (let ((env (ignore-errors
                       (genetics--json-parse
                        (buffer-substring-no-properties (point) (line-end-position))
                        :object-type 'alist :array-type 'list
                        :null-object nil :false-object nil))))
            (when (equal (alist-get 'kind env) "pipeline-run")
              (setq steps (alist-get 'data env)))))
        steps))))

(defun genetics-fastq--finished (buffer status)
  "Offer to open the VCF made by the run in BUFFER when STATUS is success."
  (when (string-prefix-p "finished" status)
    (let ((vcf (or (genetics-fastq--vcf-output
                    (genetics-fastq--result-steps buffer))
                   (buffer-local-value 'genetics-fastq--planned-vcf buffer))))
      (cond ((null vcf)
             (message "genome pipeline finished; it did not name a VCF output"))
            ((not (file-exists-p vcf))
             (message "genome pipeline finished, but %s does not exist" vcf))
            ((y-or-n-p (format "Pipeline finished; open %s? "
                               (abbreviate-file-name vcf)))
             (require 'genetics)
             (genetics-open vcf))))))

(defun genetics-fastq--start (files steps)
  "Start `genome pipeline run' on FILES asynchronously; return its buffer.
STEPS is the plan, used to find the VCF if the run does not name it."
  (let* ((argv (genetics-fastq-argv 'run files))
         (exe (executable-find (car argv))))
    (unless exe
      (genetics--error 'genetics-genome-missing
                       "genome-cli executable %S not found; install genome-cli or set `genetics-genome-executable'"
                       (car argv)))
    (let ((buf (compilation-start
                (genetics-genome-argv-string (cons exe (cdr argv)))
                #'genetics-fastq-run-mode
                (lambda (_mode) "*genetics-fastq-run*"))))
      (with-current-buffer buf
        (setq genetics-fastq--planned-vcf (genetics-fastq--vcf-output steps))
        (add-hook 'compilation-finish-functions #'genetics-fastq--finished
                  nil t))
      buf)))

(defun genetics-fastq-run (files &optional no-confirm)
  "Turn FASTQ FILES into a VCF with `genome pipeline run'.
The plan is shown first and the run starts after confirmation (skipped
when NO-CONFIRM is non-nil).  The run is asynchronous, in a compilation
buffer; when it succeeds you are offered to open the resulting VCF.
Returns the run buffer, or nil if not confirmed."
  (interactive (list (genetics--read-fastq-files)))
  (let* ((plan (genetics-fastq-plan files))
         (steps (buffer-local-value 'genetics-fastq--steps plan)))
    (when (or no-confirm
              (y-or-n-p (format "Run the pipeline on %d file(s) (it can take hours)? "
                                (length files))))
      (genetics-fastq--start files steps))))

(defun genetics-fastq-plan-run ()
  "Run the pipeline shown in this plan buffer, after confirmation."
  (interactive)
  (unless (derived-mode-p 'genetics-fastq-plan-mode)
    (user-error "Not in a FASTQ plan buffer"))
  (when (y-or-n-p "Run this pipeline (it can take hours)? ")
    (genetics-fastq--start genetics-fastq--files genetics-fastq--steps)))

(provide 'genetics-fastq)
;;; genetics-fastq.el ends here
