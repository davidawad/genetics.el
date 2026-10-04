;;; genetics-fastq-test.el --- FASTQ detection and genome pipeline commands -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(defun genetics-fastq-test--write (name)
  "Write a tiny synthetic FASTQ called NAME in a temp dir; return its path."
  (let ((file (expand-file-name name (genetics-test-temp-dir))))
    (with-temp-file file
      (insert "@SYNTH:1:FAKE:1:1101:1000:1000 1:N:0:1\nACGTNACGTA\n+\nFFFFFFFFFF\n"))
    file))

(defun genetics-fastq-test--wait (buffer)
  "Wait until the process of BUFFER has finished and its sentinel ran."
  (let ((proc (get-buffer-process buffer)) (n 0))
    (while (and proc (process-live-p proc) (< n 200))
      (accept-process-output proc 0.05)
      (cl-incf n))
    ;; let the sentinel and compilation-finish-functions run
    (accept-process-output nil 0.2)))

(ert-deftest genetics-fastq-test-detection ()
  (genetics-test-with-env
    (should (genetics-fastq-file-p "/x/sample_R1.fastq.gz"))
    (should (genetics-fastq-file-p "/x/sample.FQ"))
    (should-not (genetics-fastq-file-p (genetics-test-fixture "sample.vcf")))
    (let ((odd (make-temp-file "reads" nil ".txt")))
      (unwind-protect
          (progn
            (with-temp-file odd (insert "@r1\nACGT\n+\nFFFF\n"))
            (should (genetics-fastq-file-p odd))
            (should (eq 'fastq (genetics-detect-format odd)))
            (let ((err (should-error (genetics-parse-file odd)
                                     :type 'genetics-fastq-file)))
              (should (string-match-p "raw sequencer reads" (error-message-string err)))))
        (delete-file odd)))
    (let ((bam (make-temp-file "x" nil ".bam")))
      (unwind-protect (should (eq 'alignment (genetics-detect-format bam)))
        (delete-file bam)))))

(ert-deftest genetics-fastq-test-open-explains ()
  (genetics-test-with-env
    (let ((fq (genetics-fastq-test--write "sample_R1.fastq")))
      (should-not (genetics-open fq))
      (should-not genetics-loaded-kits)
      (let ((text (genetics-test-buffer-text "*genetics: raw reads*")))
        (should (string-match-p "raw sequencer reads, not genotypes" text))
        (should (string-match-p "open that" text))
        (should (string-match-p "\\[Show the genome-cli pipeline plan\\]" text))
        (should (string-match-p "Nothing was loaded" text))))
    (let ((bam (make-temp-file "aligned" nil ".bam")))
      (unwind-protect
          (progn
            (should-not (genetics-open bam))
            (should (string-match-p "aligned reads"
                                    (genetics-test-buffer-text "*genetics: raw reads*"))))
        (delete-file bam)))))

(ert-deftest genetics-fastq-test-argv-explain ()
  (let ((genetics-genome-executable "genome")
        (genetics-fastq-build "GRCh38")
        (genetics-fastq-output-directory "/out/")
        (genetics-fastq-extra-args '("--threads" "8")))
    (should (equal (genetics-fastq-argv 'plan '("/r/a_R1.fq.gz" "/r/a_R2.fq.gz"))
                   (list "genome" "pipeline" "plan" (expand-file-name "/r/a_R1.fq.gz")
                         (expand-file-name "/r/a_R2.fq.gz")
                         "--build" "GRCh38" "--out" (expand-file-name "/out/")
                         "--threads" "8" "--format" "json")))
    (should (equal (genetics-fastq-run-explain '("/r/a b.fq"))
                   (format "genome pipeline run %s --build GRCh38 --out %s --threads 8 --format json"
                           (shell-quote-argument (expand-file-name "/r/a b.fq"))
                           (shell-quote-argument (expand-file-name "/out/")))))
    (unless (eq system-type 'windows-nt)
      (should (equal (genetics-fastq-run-explain '("/r/a b.fq"))
                     "genome pipeline run /r/a\\ b.fq --build GRCh38 --out /out/ --threads 8 --format json")))
    (should (string-prefix-p (concat "genome pipeline plan "
                                     (shell-quote-argument (expand-file-name "/r/x.fq")))
                             (genetics-fastq-plan-explain '("/r/x.fq"))))
    (should-error (genetics-fastq-argv 'plan nil) :type 'genetics-file-error)
    (should-error (genetics-fastq-argv 'destroy '("/a.fq")) :type 'genetics-genome-error)))

(ert-deftest genetics-fastq-test-plan-buffer ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((fq (genetics-fastq-test--write "s_R1.fastq.gz"))
             (buf (genetics-fastq-plan (list fq)))
             (text (genetics-test-buffer-text buf)))
        (should (equal (list (concat "pipeline plan " fq " --build GRCh38 --format json"))
                       (funcall log)))
        (dolist (re '("1\\. fetch-reference +curl +skipped-cached"
                      "2\\. align +bwa-mem2 +planned"
                      "\\$ bwa-mem2 mem /cache/GRCh38.fa R1.fq.gz R2.fq.gz"
                      "out /out/sample.bam"
                      "Resulting VCF: /out/sample.vcf.gz"
                      "none of your data is sent"
                      "Run with: +\\$ .*pipeline run"))
          (should (string-match-p re text)))
        (with-current-buffer buf
          (should (eq major-mode 'genetics-fastq-plan-mode))
          (should (eq (key-binding "x") #'genetics-fastq-plan-run)))))))

(ert-deftest genetics-fastq-test-run-async-and-open ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((genetics-source-function #'genetics-source-native)
             (fq (genetics-fastq-test--write "s_R1.fastq.gz"))
             (asked nil))
        ;; declining runs nothing
        (cl-letf (((symbol-function 'y-or-n-p) (lambda (_p) nil)))
          (should-not (genetics-fastq-run (list fq))))
        (should-not (cl-some (lambda (l) (string-prefix-p "pipeline run" l)) (funcall log)))
        (cl-letf (((symbol-function 'y-or-n-p)
                   (lambda (p) (push p asked) t)))
          (let ((buf (genetics-fastq-run (list fq))))
            (should (equal "*genetics-fastq-run*" (buffer-name buf)))
            (with-current-buffer buf
              (should (derived-mode-p 'compilation-mode)))
            (genetics-fastq-test--wait buf)
            (should (string-match-p "\\[2/3\\] call: deepvariant"
                                    (genetics-test-buffer-text buf)))
            (should (member (concat "pipeline run " fq " --build GRCh38 --format json")
                            (funcall log)))
            ;; plan shown, run confirmed, then the VCF offered and opened
            (should (cl-some (lambda (p) (string-match-p "open .*wgs-grch38.vcf" p)) asked))
            (should (equal "wgs-grch38" (genetics-kit-name (car genetics-loaded-kits))))))))))

(ert-deftest genetics-fastq-test-run-failure-offers-nothing ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome _log
      (let* ((genetics-fastq-extra-args '("--fail"))
             (fq (genetics-fastq-test--write "s_R1.fastq.gz"))
             (asked nil))
        (cl-letf (((symbol-function 'y-or-n-p) (lambda (p) (push p asked) t)))
          (genetics-fastq-test--wait (genetics-fastq-run (list fq) t))
          (should-not asked)
          (should-not genetics-loaded-kits)
          (should (string-match-p "step_failed"
                                  (genetics-test-buffer-text "*genetics-fastq-run*"))))))))

(ert-deftest genetics-privacy-test-local-processes-only-where-expected ()
  "Only gzip (parse) and genome-cli (source, fastq) are run, all locally."
  (let ((rx "call-process\\|process-file\\|compilation-start\\|shell-command")
        (allowed '("genetics-parse.el" "genetics-source.el" "genetics-fastq.el")))
    (dolist (f (directory-files (genetics-test-root) t "\\`genetics.*\\.el\\'"))
      (unless (member (file-name-nondirectory f) allowed)
        (should-not (string-match-p rx (with-temp-buffer (insert-file-contents f)
                                                         (buffer-string))))))))

(provide 'genetics-fastq-test)
;;; genetics-fastq-test.el ends here
