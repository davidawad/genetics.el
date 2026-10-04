;;; genetics-org-test.el --- Org dynamic blocks -*- lexical-binding: t; -*-

;;; Commentary:

;; genetics-summary, genetics-hits and genetics-apoe against the
;; synthetic fixtures, natively and through the fake genome-cli.

;;; Code:

(require 'genetics-test-util)

(defun genetics-org-test--update (text &optional dir)
  "Return TEXT (an Org document) after updating all its dynamic blocks.
DIR is the document's directory (default: the fixtures directory)."
  (with-temp-buffer
    (setq default-directory (or dir (file-name-as-directory (genetics-test-fixture ""))))
    (insert text)
    (delay-mode-hooks (org-mode))
    (let ((inhibit-message t)) (org-update-all-dblocks))
    (buffer-substring-no-properties (point-min) (point-max))))

(defun genetics-org-test--block (name params)
  "Return the text written by dynamic block NAME for PARAMS."
  (with-temp-buffer
    (funcall (intern (concat "org-dblock-write:" name)) params)
    (buffer-string)))

(ert-deftest genetics-org-test-summary-file-and-reuse ()
  (genetics-test-with-env
    (let ((text (genetics-org-test--update
                 "* Kit\n#+BEGIN: genetics-summary :file \"23andme-sample.txt\"\n#+END:\n")))
      (dolist (re '("^- Format :: 23andMe$" "^- Assay :: array, every assayed site listed$"
                    "^- Build :: GRCh37$" "^- Records :: 33$"
                    "^- No-call rate :: 18.18% (6 no-calls)$" "^- Inferred sex :: female"
                    "^- Caveats ::$" "^  - .*no liftover" "Palindromic SNPs"
                    "Informational only, not medical advice" "^#\\+END:$"))
        (should (string-match-p re text)))
      ;; registered, but no buffer was shown
      (should (= 1 (length genetics-loaded-kits)))
      (should-not (get-buffer "*genetics: 23andme-sample*")))
    ;; a second block on the same file reuses the loaded kit
    (let ((kit (car genetics-loaded-kits)))
      (genetics-org-test--block "genetics-summary"
                                (list :file (genetics-test-fixture "23andme-sample.txt")))
      (should (equal (list kit) genetics-loaded-kits))
      (should (string-match-p "Reuse the kit \"23andme-sample\""
                              (genetics-org-summary-explain
                               (list :file (genetics-test-fixture "23andme-sample.txt"))))))))

(ert-deftest genetics-org-test-hits-table-and-filters ()
  (genetics-test-with-env
    (genetics-test-load "23andme-sample.txt")
    (let ((text (genetics-org-test--block "genetics-hits" '(:kit "23andme-sample"))))
      (should (string-match-p "^| rsid +| gene +| genotype +| risk-allele copies +| call source +| magnitude +| effect +| source +|$" text))
      (should (= 8 (cl-count-if (lambda (l) (string-match-p "^| rs[0-9]+ " l))
                                (split-string text "\n"))))
      (should (string-match-p "^| rs1801133 +| MTHFR +| AG +| 1 (A) +| observed +| +2 | MTHFR C677T (A222V)\\. +| \\[\\[https://www.snpedia.com/index.php/Rs1801133\\]\\[source\\]\\] |$" text))
      (should (string-match-p "never flipped automatically" text))
      (should (string-match-p "Informational only, not medical advice" text)))
    (let ((text (genetics-org-test--block
                 "genetics-hits" '(:kit "23andme-sample" :min-magnitude 3 :genes (apoe "F5")))))
      (should (equal '("rs429358" "rs6025")
                     (sort (delq nil (mapcar (lambda (l) (and (string-match "^| \\(rs[0-9]+\\) " l)
                                                               (match-string 1 l)))
                                             (split-string text "\n")))
                           #'string<)))
      (should (string-match-p "Filter: magnitude >= 3; genes APOE, F5\\." text)))
    ;; :effect-width adds an Org width cookie row for the effect column
    (should (string-match-p "^| +| +| +| +| +| +| <40> +| +|$"
                            (genetics-org-test--block
                             "genetics-hits" '(:kit "23andme-sample" :effect-width 40))))
    ;; a string magnitude is accepted
    (should (string-match-p "Filter: magnitude >= 4\\."
                            (genetics-org-test--block
                             "genetics-hits" '(:kit "23andme-sample" :min-magnitude "4"))))
    (should (string-match-p "passes the filter (8 before filtering)"
                            (genetics-org-test--block "genetics-hits"
                                                      '(:kit "23andme-sample" :genes "NOPE"))))
    (should (string-match-p "magnitude >= 2 in MTHFR"
                            (genetics-org-hits-explain '(:min-magnitude 2 :genes "mthfr"))))))

(ert-deftest genetics-org-test-strand-flip-note ()
  (genetics-test-with-env
    (genetics-test-kit-from-rows '(("rs1801133" "1" 11856378 "TT")))
    (let ((text (genetics-org-test--block "genetics-hits" nil)))
      (should (string-match-p "| rs1801133 +| MTHFR +| TT +| 0 (A) +| observed " text))
      (should (string-match-p "^- rs1801133 :: Possible strand flip" text)))
    (should (string-match-p "could not be assessed"
                            (genetics-org-test--block "genetics-apoe" nil)))))

(ert-deftest genetics-org-test-wgs-inferred-native ()
  (genetics-test-with-env
    (genetics-test-load "wgs-grch38.vcf" :ref-calls 'absent-means-ref)
    (let ((hits (genetics-org-test--block "genetics-hits" nil))
          (apoe (genetics-org-test--block "genetics-apoe" nil))
          (summary (genetics-org-test--block "genetics-summary" nil)))
      (should (string-match-p "| rs429358 +| APOE +| TT +| 0 (C) +| inferred ref +|" hits))
      (should (string-match-p "Call source \"inferred ref\": Inferred, not observed" hits))
      (should (string-match-p "\\*APOE: e2/e3\\* (status: ok)" apoe))
      (should (string-match-p "rs429358 TT (inferred ref), rs7412 CT" apoe))
      (should (string-match-p "inferred homozygous reference" apoe))
      (should (string-match-p "+ strand" apoe))
      (should (string-match-p "Informational only, not medical advice" apoe))
      (should (string-match-p "^- Assay :: wgs, variant sites only (absent = reference, inferred)$" summary))
      (should (string-match-p "^- Build :: GRCh38$" summary)))))

(ert-deftest genetics-org-test-fake-genome ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((file (genetics-test-fixture "wgs-grch38.vcf"))
             (explain (genetics-org-summary-explain (list :file file))))
        ;; the explain twin runs nothing
        (should (string-match-p "Open .*wgs-grch38.vcf with genome-cli" explain))
        (should (string-match-p (regexp-quote (genetics-source-genome-cli-explain file)) explain))
        (should-not (funcall log))
        (let ((text (genetics-org-test--update
                     (concat "#+BEGIN: genetics-summary :file \"wgs-grch38.vcf\"\n#+END:\n"
                             "#+BEGIN: genetics-hits :file \"wgs-grch38.vcf\" :genes (APOE MTHFR)\n#+END:\n"
                             "#+BEGIN: genetics-apoe :kit \"wgs-grch38\"\n#+END:\n"))))
          ;; imported once; the later blocks reuse the registered kit
          (should (= 1 (cl-count-if (lambda (l) (string-prefix-p "import " l)) (funcall log))))
          (should (= 1 (length genetics-loaded-kits)))
          (dolist (re '("^- Source :: genome-cli (kit k1)$" "^- Sample :: SYNTH01$"
                        "^- Records :: 27$" "^- Inferred sex :: male" "synthetic fixture"
                        "| rs429358 +| APOE +| TT +| 0 (C) +| inferred ref +|"
                        "| rs7412 +| APOE +| CT +| 1 (T) +| observed +|"
                        "| rs1801133 +| MTHFR +| GA +| 1 (A) +| observed +|"
                        "\\*APOE: e2/e3\\*" "rs429358 TT (inferred ref)"))
            (should (string-match-p re text)))
          (should (= 3 (cl-count-if (lambda (l) (string-match-p "Informational only, not medical advice" l))
                                    (split-string text "\n")))))))))

(ert-deftest genetics-org-test-failures-are-comments ()
  (genetics-test-with-env
    ;; no kit loaded at all
    (let ((text (genetics-org-test--block "genetics-summary" nil)))
      (should (string-match-p "\\`# genetics-summary failed: No genetics kit available: No kit loaded" text))
      (should (string-match-p "^# What to do: Load the kit with M-x genetics-open" text))
      (should (string-match-p "Informational only, not medical advice" text)))
    ;; missing file
    (let ((text (genetics-org-test--block "genetics-hits" '(:file "/nonexistent/kit.txt"))))
      (should (string-match-p "\\`# genetics-hits failed: .*Cannot read /nonexistent/kit.txt" text))
      (should (string-match-p "^# What to do: Check that :file" text)))
    ;; unknown kit name; the explain twin predicts the failure
    (should (string-match-p "^# genetics-apoe failed: .*No loaded kit named \"nope\""
                            (genetics-org-test--block "genetics-apoe" '(:kit "nope"))))
    (should (string-match-p "not loaded; the block would fail"
                            (genetics-org-apoe-explain '(:kit "nope"))))
    ;; genome-cli failure
    (genetics-test-with-fake-genome _log
      (let ((text (genetics-org-test--block
                   "genetics-summary" (list :file (genetics-test-fixture "missing-kit.vcf")))))
        (should (string-match-p "^# genetics-summary failed" text)))
      (let* ((file (make-temp-file "genetics-missing-" nil ".vcf"))
             (text (unwind-protect
                       (genetics-org-test--block "genetics-summary" (list :file file))
                     (delete-file file))))
        (should (string-match-p "^# genetics-summary failed: genome-cli failed" text))
        (should (string-match-p "^# What to do: Run the command shown by M-x genetics-source-genome-cli-explain" text))))
    ;; several kits and no :kit
    (genetics-test-load "23andme-sample.txt")
    (genetics-test-load "myheritage-sample.csv")
    (should (string-match-p "2 kits are loaded"
                            (genetics-org-test--block "genetics-apoe" nil)))
    ;; whatever fails, the document survives: every non-comment line is intact
    (let ((text (genetics-org-test--update
                 "* A\n#+BEGIN: genetics-apoe :kit \"nope\"\n#+END:\n* B\ntext\n")))
      (should (string-match-p "\\`\\* A\n#\\+BEGIN: genetics-apoe :kit \"nope\"\n# genetics-apoe failed" text))
      (should (string-suffix-p "#+END:\n* B\ntext\n" text)))))

(ert-deftest genetics-org-test-explain-block-and-insert ()
  (genetics-test-with-env
    (genetics-test-load "23andme-sample.txt")
    (with-temp-buffer
      (insert "#+BEGIN: genetics-hits :kit \"23andme-sample\" :min-magnitude 3\nold\n#+END:\n")
      (delay-mode-hooks (org-mode))
      (goto-char (point-min))
      (forward-line 1)
      (should (string-match-p "\\`genetics-hits: Use the loaded kit \"23andme-sample\". .*magnitude >= 3"
                              (genetics-org-explain-block)))
      ;; explaining does not touch the block
      (should (string-match-p "^old$" (buffer-string))))
    (with-temp-buffer
      (delay-mode-hooks (org-mode))
      (let ((inhibit-message t))
        (genetics-org-insert-block "genetics-apoe" "23andme-sample"))
      (should (string-match-p "\\`#\\+BEGIN: genetics-apoe :kit \"23andme-sample\"\n\\*APOE: e2/e4\\*"
                              (buffer-string))))))

(provide 'genetics-org-test)
;;; genetics-org-test.el ends here
