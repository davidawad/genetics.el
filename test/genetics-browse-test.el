;;; genetics-browse-test.el --- browser filters, lookup and summary -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(defun genetics-browse-test--ids (kit filters &optional limit)
  "Return rsids of KIT matching FILTERS (with LIMIT)."
  (mapcar #'genetics-snp-rsid (car (genetics-browse-rows kit filters limit))))

(ert-deftest genetics-browse-test-filters ()
  (genetics-test-with-env
    (let ((kit (genetics-test-load "23andme-sample.txt")))
      (should (= 33 (length (genetics-browse-test--ids kit nil))))
      (should (= 10 (length (genetics-browse-test--ids kit '(:chrom "X")))))
      (should (equal '("rs3094315" "rs3131972" "rs12124819")
                     (genetics-browse-test--ids kit '(:chrom "1" :range (752566 . 776546)))))
      (should (= 10 (length (genetics-browse-test--ids kit '(:rsid "^rs50")))))
      (should (= 1 (length (genetics-browse-test--ids kit '(:rsid "RS4477212")))))
      ;; genotype filter ignores allele order (key of "GA" == key of "AG")
      (let ((n (length (genetics-browse-test--ids
                        kit (list :genotype (genetics-genotype-key "GA"))))))
        (should (= 5 n)))
      (should (= 6 (length (genetics-browse-test--ids kit '(:zygosity no-call)))))
      (should (= 10 (length (genetics-browse-test--ids kit '(:zygosity heterozygous)))))
      (should (= 15 (length (genetics-browse-test--ids kit '(:zygosity homozygous)))))
      (should (= 8 (length (genetics-browse-test--ids kit '(:annotated t)))))
      (should (= 2 (length (genetics-browse-test--ids kit '(:annotated t :chrom "19" :genotype "C/T")))))
      ;; combined filters
      (should (equal '("rs5001" "rs5003" "rs5005" "rs5008")
                     (genetics-browse-test--ids kit '(:chrom "X" :zygosity heterozygous)))))))

(ert-deftest genetics-browse-test-limit-truncates ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (res (genetics-browse-rows kit nil 5)))
      (should (= 5 (length (car res))))
      (should (cdr res))
      (should-not (cdr (genetics-browse-rows kit nil 33)))
      (should-not (cdr (genetics-browse-rows kit nil nil))))))

(ert-deftest genetics-browse-test-buffer-and-commands ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (buf (genetics-browse kit)))
      (with-current-buffer buf
        (should (derived-mode-p 'genetics-browse-mode))
        (should (eq kit genetics--buffer-kit))
        (should (= 33 (count-lines (point-min) (point-max))))
        (should (string-match-p "Filters: none" header-line-format))
        (genetics-browse-filter-chromosome "chrX")
        (should (= 10 (count-lines (point-min) (point-max))))
        (should (string-match-p "chr=X" header-line-format))
        (genetics-browse-filter-heterozygous)
        (should (= 4 (count-lines (point-min) (point-max))))
        (genetics-browse-filter-no-calls)
        (should (= 0 (count-lines (point-min) (point-max))))
        (genetics-browse-clear-filters)
        (should (= 33 (count-lines (point-min) (point-max))))
        (genetics-browse-filter-annotated)
        (should (= 8 (count-lines (point-min) (point-max))))
        (genetics-browse-filter-annotated)
        (should (= 33 (count-lines (point-min) (point-max))))
        (genetics-browse-filter-range 100 250)
        (should (= 4 (count-lines (point-min) (point-max))))
        (genetics-browse-clear-filters)
        (genetics-browse-filter-rsid "^rs44")
        (should (= 1 (count-lines (point-min) (point-max))))
        (genetics-browse-clear-filters)
        (genetics-browse-filter-genotype "ga")
        (should (= 5 (count-lines (point-min) (point-max))))
        (should (string-match-p "genotype=A/G" header-line-format))
        (should-error (genetics-browse-filter-range 10 5) :type 'user-error)
        (should-error (genetics-browse-filter-rsid "[") :type 'invalid-regexp)
        (genetics-browse-clear-filters)
        ;; gene column comes from annotations
        (goto-char (point-min))
        (search-forward "rs1801133")
        (should (string-match-p "MTHFR" (buffer-substring (line-beginning-position) (line-end-position))))))))

(ert-deftest genetics-browse-test-limit-in-buffer ()
  (genetics-test-with-env
    (let* ((genetics-browse-limit 5)
           (kit (genetics-test-load "23andme-sample.txt"))
           (buf (genetics-browse kit)))
      (with-current-buffer buf
        (should (= 5 (count-lines (point-min) (point-max))))
        (should genetics-browse--truncated)
        (should (string-match-p "truncated" header-line-format))))))

(ert-deftest genetics-browse-test-ret-opens-lookup ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (buf (genetics-browse kit)))
      (with-current-buffer buf
        (goto-char (point-min))
        (search-forward "rs7412")
        (genetics-browse-lookup))
      (should (get-buffer "*genetics-lookup: rs7412*")))))

(ert-deftest genetics-browse-test-filter-commands-need-browser ()
  (with-temp-buffer
    (should-error (genetics-browse-clear-filters) :type 'user-error)))

(ert-deftest genetics-lookup-test-buffer-contents ()
  (genetics-test-with-env
    (genetics-test-load "23andme-sample.txt")
    (genetics-test-load "ancestry-sample.txt")
    (let ((buf (genetics-lookup " RS1801133 ")))
      (should (equal "*genetics-lookup: rs1801133*" (buffer-name buf)))
      (let ((text (genetics-test-buffer-text buf)))
        (should (string-match-p "23andme-sample +AG" text))
        (should (string-match-p "ancestry-sample +GG" text))
        (should (string-match-p "MTHFR" text))
        (should (string-match-p "Risk allele: +A" text))
        (should (string-match-p "One copy of 677T" text))
        (should (string-match-p "No copies of 677T" text))
        (should (string-match-p "Strand:" text))
        (should (string-match-p "snpedia.com" text))
        (should-not (string-match-p "SNPedia\n" text))))
    (let ((text (genetics-test-buffer-text (genetics-lookup "rs429358"))))
      (should (string-match-p "APOE haplotype in 23andme-sample: e2/e4" text))
      (should (string-match-p "APOE haplotype in ancestry-sample: e3/e4" text)))
    (let ((text (genetics-test-buffer-text (genetics-lookup "rs11240777"))))
      (should (string-match-p "No annotation" text))
      (should (string-match-p "ancestry-sample +not present" text)))))

(ert-deftest genetics-lookup-test-strand-flip-shown ()
  (genetics-test-with-env
    (genetics-test-kit-from-rows '(("rs1801133" "1" 11856378 "TT")))
    (let ((text (genetics-test-buffer-text (genetics-lookup "rs1801133"))))
      (should (string-match-p "Flag: +strand-flip" text))
      (should (string-match-p "Possible strand flip" text)))))

(ert-deftest genetics-lookup-test-requires-kit ()
  (genetics-test-with-env
    (should-error (genetics-lookup "rs1") :type 'genetics-no-kit)))

(ert-deftest genetics-summary-test-open-and-summary-buffer ()
  (genetics-test-with-env
    (let ((kit (genetics-open (genetics-test-fixture "23andme-sample.txt"))))
      (should (memq kit genetics-loaded-kits))
      (let ((text (genetics-test-buffer-text "*genetics: 23andme-sample*")))
        (dolist (re '("Format: +23andMe" "Build: +GRCh37" "Chip: +unknown"
                      "SNPs: +33" "No-calls: +6 (18\\.18%)"
                      "10 heterozygous, 15 homozygous, 2 hemizygous"
                      "Sex: +female" "Per-chromosome counts" "X +10"
                      "forward" "no liftover" "not medical advice"
                      "\\[Browse\\]" "\\[Report\\]" "\\[Lookup\\]"))
          (should (string-match-p re text))))
      (with-current-buffer "*genetics: 23andme-sample*"
        (should (derived-mode-p 'genetics-summary-mode))
        (should (eq kit genetics--buffer-kit))
        (should (eq 'genetics-summary-browse (lookup-key genetics-summary-mode-map "b")))
        (genetics-summary-browse)
        (should (get-buffer "*genetics-browse: 23andme-sample*"))
        (genetics-summary-report)
        (should (get-buffer "*genetics-report: 23andme-sample*")))
      ;; reopening the same file replaces the kit; a copy gets a unique name
      (genetics-open (genetics-test-fixture "23andme-sample.txt"))
      (should (= 1 (length genetics-loaded-kits)))
      (genetics-close "23andme-sample")
      (should (null genetics-loaded-kits))
      (should-not (get-buffer "*genetics: 23andme-sample*")))))

(ert-deftest genetics-summary-test-build-mismatch-warning-and-names ()
  (genetics-test-with-env
    (let ((warned nil) (f38 (genetics-test-write-23andme
                             '(("rs1" "1" 10 "AA"))
                             "# reference human assembly build 38\n")))
      (unwind-protect
          (cl-letf (((symbol-function 'display-warning)
                     (lambda (_type msg &rest _) (setq warned msg))))
            (genetics-open (genetics-test-fixture "23andme-sample.txt"))
            (should-not warned)
            (genetics-open f38)
            (should (string-match-p "different genome builds" warned)))
        (delete-file f38)))))

(ert-deftest genetics-summary-test-vcf-lazy-summary ()
  (genetics-test-with-env
    (let ((genetics-vcf-eager-limit 10))
      (genetics-open (genetics-test-fixture "sample.vcf"))
      (let ((text (genetics-test-buffer-text "*genetics: sample*")))
        (should (string-match-p "offset-indexed" text))
        (should (string-match-p "REF/ALT" text))))))

(ert-deftest genetics-summary-test-data-directory-default ()
  (should (equal "~/Documents/Genetics/"
                 (default-value 'genetics-data-directory))))

(provide 'genetics-browse-test)
;;; genetics-browse-test.el ends here
