;;; genetics-report-test.el --- report, compare and export -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)
(require 'json)

(ert-deftest genetics-report-test-contents ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (text (genetics-report-string kit)))
      (should (string-prefix-p "#+TITLE: Genetics report: 23andme-sample" text))
      (dolist (re '("^\\* Kit$" "Format: 23andMe" "Build: GRCh37" "Records: 33"
                    "^\\* Strand and build caveats$" "(forward) strand of GRCh37"
                    "no liftover" "never flipped automatically"
                    "^\\* Annotated findings$"
                    "| rsid | gene | genotype | risk allele | copies | interpretation | strand note | source |"
                    "^\\* APOE haplotype$" "Diplotype: e2/e4" "Status: ambiguous" "e1/e3"
                    "^\\* Disclaimer$" "not medical advice"))
        (should (string-match-p re text)))
      ;; one table row per annotated hit (8), with a source link
      (should (= 8 (with-temp-buffer
                     (insert text)
                     (goto-char (point-min))
                     (let ((n 0)) (while (re-search-forward "^| rs[0-9]+ |" nil t) (cl-incf n)) n))))
      (should (string-match-p "| rs1801133 | MTHFR | AG | A | 1 | " text))
      (should (string-match-p "\\[\\[https://www.snpedia.com/index.php/Rs1801133\\]\\[source\\]\\]" text)))))

(ert-deftest genetics-report-test-strand-flip-and-empty ()
  (genetics-test-with-env
    (let ((text (genetics-report-string
                 (genetics-test-kit-from-rows '(("rs1801133" "1" 11856378 "TT"))))))
      (should (string-match-p "Possible strand flip" text))
      (should (string-match-p "| rs1801133 | MTHFR | TT | A | 0 |" text))
      (should (string-match-p "could not be assessed" text)))
    (let ((text (genetics-report-string
                 (genetics-test-kit-from-rows '(("rs1" "1" 1 "AA"))))))
      (should (string-match-p "No annotated SNPs" text)))))

(ert-deftest genetics-report-test-buffer-and-file ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (file (make-temp-file "genetics-test-" nil ".org")))
      (unwind-protect
          (let ((buf (genetics-report kit file)))
            (with-current-buffer buf (should (derived-mode-p 'org-mode)))
            (should (equal (genetics-test-buffer-text buf)
                           (with-temp-buffer (insert-file-contents file) (buffer-string))))
            (should (equal (genetics-report-write kit file) file)))
        (delete-file file)))))

(ert-deftest genetics-compare-test-fixtures ()
  (genetics-test-with-env
    (let* ((a (genetics-test-load "23andme-sample.txt"))
           (b (genetics-test-load "myheritage-sample.csv"))
           (r (genetics-compare-kits a b)))
      (should (= 8 (plist-get r :overlap)))
      (should (= 7 (plist-get r :compared)))       ; rs3131972 is a no-call in both
      (should (= 5 (plist-get r :concordant)))
      (should (equal '("rs429358" "rs7412") (sort (mapcar #'car (plist-get r :discordant)) #'string<)))
      (should (< (abs (- 71.43 (plist-get r :concordance))) 0.01))
      (should-not (plist-get r :build-mismatch))
      (let ((text (genetics-test-buffer-text (genetics-compare a "myheritage-sample"))))
        (should (string-match-p "Overlapping ids: +8" text))
        (should (string-match-p "Concordance: +71\\.43%" text))
        (should (string-match-p "rs429358 +CT +CC" text))
        (should-not (string-match-p "WARNING" text))))))

(ert-deftest genetics-compare-test-order-hemizygous-and-complement ()
  (genetics-test-with-env
    (let* ((a (genetics-test-kit-from-rows
               '(("rs1" "1" 1 "AG") ("rs2" "1" 2 "AA") ("rs3" "X" 3 "A") ("rs4" "1" 4 "AT")
                 ("rs5" "1" 5 "CG") ("rs6" "X" 6 "A") ("rs7" "1" 7 "--") ("rs8" "1" 8 "AA"))))
           (b (genetics-test-kit-from-rows
               '(("rs1" "1" 1 "GA") ("rs2" "1" 2 "TT") ("rs3" "X" 3 "AA") ("rs4" "1" 4 "TA")
                 ("rs5" "1" 5 "GC") ("rs6" "X" 6 "AG") ("rs7" "1" 7 "AA") ("rs9" "1" 9 "AA"))))
           (r (genetics-compare-kits a b)))
      (should (= 7 (plist-get r :overlap)))
      (should (= 6 (plist-get r :compared)))        ; rs7 no-call in A
      ;; AG==GA, A==AA (hemizygous), AT==TA, CG==GC
      (should (= 4 (plist-get r :concordant)))
      ;; AA vs TT is a complement-strand difference, reported separately
      (should (equal '(("rs2" "AA" "TT")) (plist-get r :complement)))
      ;; hemizygous A vs heterozygous AG is a real discordance
      (should (equal '(("rs6" "A" "AG")) (plist-get r :discordant)))
      (should (< (abs (- (plist-get r :concordance) 66.667)) 0.01)))))

(ert-deftest genetics-compare-test-build-mismatch-and-empty ()
  (genetics-test-with-env
    (let* ((a (genetics-test-kit-from-rows '(("rs1" "1" 1 "AA"))))
           (b (genetics-test-kit-from-rows '(("rs1" "1" 5 "AA"))
                                           "# reference human assembly build 38\n"))
           (c (genetics-test-kit-from-rows '(("rs77" "1" 5 "AA")))))
      (should (plist-get (genetics-compare-kits a b) :build-mismatch))
      (should (string-match-p "different genome builds"
                              (genetics-test-buffer-text (genetics-compare a b))))
      (let ((r (genetics-compare-kits a c)))
        (should (= 0 (plist-get r :overlap)))
        (should (null (plist-get r :concordance))))
      (should-error (genetics-compare a "nope") :type 'genetics-no-kit))))

(ert-deftest genetics-export-test-csv-round-trip ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (file (make-temp-file "genetics-test-" nil ".csv")))
      (unwind-protect
          (progn
            (should (= 33 (genetics-export-csv file kit)))
            (let ((text (with-temp-buffer (insert-file-contents file) (buffer-string))))
              (should (string-prefix-p "RSID,CHROMOSOME,POSITION,RESULT,ZYGOSITY,GENE,REF,ALT\n" text))
              (should (string-match-p "\"rs1801133\",\"1\",\"11856378\",\"AG\",\"heterozygous\",\"MTHFR\"" text)))
            ;; the exported CSV is itself a readable genotype file
            (let ((back (genetics-parse-file file)))
              (should (eq 'ftdna (genetics-kit-format back)))
              (should (equal (mapcar (lambda (s) (list (genetics-snp-rsid s) (genetics-snp-chrom s)
                                                       (genetics-snp-pos s) (genetics-snp-genotype s)))
                                     (genetics-test-snps kit))
                             (mapcar (lambda (s) (list (genetics-snp-rsid s) (genetics-snp-chrom s)
                                                       (genetics-snp-pos s) (genetics-snp-genotype s)))
                                     (genetics-test-snps back)))))
            (should (= 10 (genetics-export-csv file kit '(:chrom "X")))))
        (delete-file file)))))

(ert-deftest genetics-export-test-json-round-trip ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "sample.vcf"))
           (file (make-temp-file "genetics-test-" nil ".json")))
      (unwind-protect
          (progn
            (should (= 13 (genetics-export-json file kit)))
            (let* ((data (with-temp-buffer
                           (insert-file-contents file)
                           (json-parse-buffer :object-type 'alist :array-type 'list :null-object nil)))
                   (row (cl-find "rs9999" data :key (lambda (r) (alist-get 'rsid r)) :test #'equal)))
              (should (= 13 (length data)))
              (should (equal "GT" (alist-get 'genotype row)))
              (should (equal 900100 (alist-get 'position row)))
              (should (equal "heterozygous" (alist-get 'zygosity row)))
              (should (equal "G,T" (alist-get 'alt row)))
              (should (null (alist-get 'gene row)))
              (should (equal "APOE" (alist-get 'gene (cl-find "rs429358" data
                                                              :key (lambda (r) (alist-get 'rsid r))
                                                              :test #'equal))))))
        (delete-file file)))))

(ert-deftest genetics-export-test-uses-browser-filters ()
  (genetics-test-with-env
    (let* ((kit (genetics-test-load "23andme-sample.txt"))
           (buf (genetics-browse kit))
           (csv (make-temp-file "genetics-test-" nil ".csv"))
           (json (make-temp-file "genetics-test-" nil ".json"))
           (genetics-browse-limit 3))
      (unwind-protect
          (with-current-buffer buf
            (genetics-browse-filter-chromosome "X")
            (genetics-browse-filter-heterozygous)
            ;; exports ignore the display limit
            (should (= 4 (genetics-export-csv csv)))
            (should (= 4 (genetics-browse-export-json json)))
            (genetics-browse-clear-filters)
            (should (= 33 (genetics-browse-export-csv csv))))
        (delete-file csv) (delete-file json)))))

(provide 'genetics-report-test)
;;; genetics-report-test.el ends here
