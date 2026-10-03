;;; genetics-annotate-test.el --- annotations, risk alleles, APOE -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(defun genetics-annotate-test--curated ()
  "Return the curated annotation file path."
  (expand-file-name "annotations/genetics-curated.json" (genetics-test-root)))

(ert-deftest genetics-annotate-test-default-files-include-curated ()
  (should (file-exists-p (car (default-value 'genetics-annotation-files))))
  (should (string-suffix-p "annotations/genetics-curated.json"
                           (car (default-value 'genetics-annotation-files)))))

(ert-deftest genetics-annotate-test-load-curated-json ()
  (let* ((anns (genetics-load-annotation-file (genetics-annotate-test--curated)))
         (ids (mapcar #'genetics-annotation-rsid anns)))
    (dolist (id '("rs429358" "rs7412" "rs1801133" "rs1801131" "rs6025"
                  "rs1800562" "rs4988235" "rs4244285"))
      (should (member id ids)))
    (let ((mthfr (cl-find "rs1801133" anns :key #'genetics-annotation-rsid :test #'equal)))
      (should (equal "MTHFR" (genetics-annotation-gene mthfr)))
      (should (equal "A" (genetics-annotation-risk-allele mthfr)))
      (should (string-match-p "gene strand" (genetics-annotation-strand mthfr)))
      (should (string-match-p "snpedia" (genetics-annotation-url mthfr)))
      (should (string-match-p "ncbi.nlm.nih.gov/snp/rs1801133" (genetics-annotation-notes mthfr)))
      (should (assoc "A/G" (genetics-annotation-genotypes mthfr))))
    (dolist (a anns)
      (should (genetics-annotation-url a))
      (should (string-match-p "not medical advice\\|Informational only" (genetics-annotation-notes a))))
    (let ((f5 (cl-find "rs6025" anns :key #'genetics-annotation-rsid :test #'equal)))
      (should (equal "T" (genetics-annotation-risk-allele f5)))
      (should (string-match-p "gene strand" (genetics-annotation-strand f5))))
    (should (equal "A" (genetics-annotation-risk-allele
                        (cl-find "rs1800562" anns :key #'genetics-annotation-rsid :test #'equal))))
    (should (equal "A" (genetics-annotation-risk-allele
                        (cl-find "rs4244285" anns :key #'genetics-annotation-rsid :test #'equal))))))

(ert-deftest genetics-annotate-test-load-org ()
  (let ((anns (genetics-load-annotation-file
               (expand-file-name "annotations/genetics-example.org" (genetics-test-root)))))
    (should (= 2 (length anns)))
    (let ((mc1r (car anns)) (herc (cadr anns)))
      (should (equal "rs1805007" (genetics-annotation-rsid mc1r)))
      (should (equal "MC1R" (genetics-annotation-gene mc1r)))
      (should (equal "T" (genetics-annotation-risk-allele mc1r)))
      (should (= 2 (genetics-annotation-magnitude mc1r)))
      (should (equal "One R151C allele." (cdr (assoc "C/T" (genetics-annotation-genotypes mc1r)))))
      (should (string-match-p "Example entry" (genetics-annotation-notes mc1r)))
      ;; rsid taken from the heading when there is no RSID property
      (should (equal "rs12913832" (genetics-annotation-rsid herc)))
      (should (string-match-p "blue eyes" (genetics-annotation-notes herc))))))

(ert-deftest genetics-annotate-test-user-files-override-and-reload ()
  (genetics-test-with-env
    (let ((user (make-temp-file "genetics-test-" nil ".json")))
      (unwind-protect
          (let ((genetics-annotation-files (list (genetics-annotate-test--curated) user)))
            (with-temp-file user
              (insert "[{\"rsid\":\"rs1801133\",\"gene\":\"MYGENE\",\"risk_allele\":\"G\"},{\"rsid\":\"RS99\",\"gene\":\"X\"}]"))
            (should (equal "MYGENE" (genetics-annotation-gene (genetics-annotation "rs1801133"))))
            (should (equal "X" (genetics-annotation-gene (genetics-annotation "rs99"))))
            (should (genetics-annotation "rs429358"))
            ;; editing the file is picked up without a reload command
            (sleep-for 0.01)
            (with-temp-file user (insert "[{\"rsid\":\"rs1801133\",\"gene\":\"CHANGED12\"}]"))
            (should (equal "CHANGED12" (genetics-annotation-gene (genetics-annotation "rs1801133"))))
            (should (null (genetics-annotation "rs99"))))
        (delete-file user)))))

(ert-deftest genetics-annotate-test-errors ()
  (let ((f (make-temp-file "genetics-test-" nil ".json")))
    (unwind-protect
        (progn
          (with-temp-file f (insert "[{\"gene\":\"X\"}]"))
          (let ((err (should-error (genetics-load-annotation-file f) :type 'genetics-annotation-error)))
            (should (string-match-p "no rsid" (cadr err))))
          (with-temp-file f (insert "not json"))
          (should-error (genetics-load-annotation-file f) :type 'genetics-annotation-error)
          (with-temp-file f (insert "{\"a\":1}"))
          (should-error (genetics-load-annotation-file f) :type 'genetics-annotation-error))
      (delete-file f))
    (should-error (genetics-load-annotation-file "/nonexistent/ann.json")
                  :type 'genetics-annotation-error)
    (let ((txt (make-temp-file "genetics-test-" nil ".txt")))
      (unwind-protect
          (should-error (genetics-load-annotation-file txt) :type 'genetics-annotation-error)
        (delete-file txt)))))

(ert-deftest genetics-annotate-test-risk-copies ()
  (should (equal 0 (plist-get (genetics-risk-assess "GG" "A" "G") :copies)))
  (should (equal 1 (plist-get (genetics-risk-assess "AG" "A" "G") :copies)))
  (should (equal 1 (plist-get (genetics-risk-assess "GA" "A" "G") :copies)))
  (should (equal 2 (plist-get (genetics-risk-assess "AA" "A" "G") :copies)))
  (should (equal 1 (plist-get (genetics-risk-assess "A" "A" "G") :copies)))
  (should (null (plist-get (genetics-risk-assess "AG" "A" "G") :flag)))
  (should (null (plist-get (genetics-risk-assess "AA" "A") :flag)))
  (let ((r (genetics-risk-assess "--" "A" "G")))
    (should (null (plist-get r :copies)))
    (should (eq 'no-call (plist-get r :flag)))))

(ert-deftest genetics-annotate-test-strand-flip-flag ()
  (let ((r (genetics-risk-assess "TT" "A" "G")))
    (should (eq 'strand-flip (plist-get r :flag)))
    (should (= 0 (plist-get r :copies)))      ; never silently flipped
    (should (= 2 (plist-get r :flipped-copies))))
  (let ((r (genetics-risk-assess "TC" "A" "G")))
    (should (eq 'strand-flip (plist-get r :flag)))
    (should (= 1 (plist-get r :flipped-copies))))
  (should (eq 'strand-flip (plist-get (genetics-risk-assess "TT" "A") :flag)))
  ;; a genuine non-risk genotype is not a flip
  (should (null (plist-get (genetics-risk-assess "GG" "A") :flag)))
  (should (eq 'unexpected (plist-get (genetics-risk-assess "AT" "A" "G") :flag))))

(ert-deftest genetics-annotate-test-palindromic-flag ()
  (should (eq 'ambiguous (plist-get (genetics-risk-assess "AT" "A" "T") :flag)))
  (should (eq 'ambiguous (plist-get (genetics-risk-assess "AA" "A" "T") :flag)))
  (should (eq 'ambiguous (plist-get (genetics-risk-assess "CG" "C" "G") :flag)))
  ;; no other allele known: only a het containing both complements is ambiguous
  (should (eq 'ambiguous (plist-get (genetics-risk-assess "AT" "A") :flag)))
  (should (= 1 (plist-get (genetics-risk-assess "AT" "A" "T") :copies))))

(ert-deftest genetics-annotate-test-assess-interpretation ()
  (genetics-test-with-env
    (let ((ann (genetics-annotation "rs1801133")))
      (let ((a (genetics-assess ann "AG")))
        (should (= 1 (plist-get a :copies)))
        (should (string-match-p "One copy of 677T" (plist-get a :interpretation))))
      (should (string-match-p "Two copies" (plist-get (genetics-assess ann "AA") :interpretation)))
      (should (string-match-p "No copies" (plist-get (genetics-assess ann "GG") :interpretation)))
      (let ((a (genetics-assess ann "TT")))
        (should (eq 'strand-flip (plist-get a :flag)))
        (should (string-match-p "Possible strand flip" (plist-get a :interpretation)))
        (should-not (string-match-p "Two copies of 677T" (plist-get a :interpretation))))
      (should (string-match-p "No call" (plist-get (genetics-assess ann "--") :interpretation))))
    ;; generic text when no per-genotype entry exists
    (let ((ann (genetics-annotation-create :rsid "rs1" :risk-allele "A" :effect "Some effect")))
      (should (string-match-p "One copy of the A allele"
                              (plist-get (genetics-assess ann "AG") :interpretation))))))

(ert-deftest genetics-annotate-test-apoe-all-genotype-combinations ()
  (let ((cases '(("TT" "CC" "e3/e3" ok) ("TT" "CT" "e2/e3" ok) ("TT" "TT" "e2/e2" ok)
                 ("CT" "CC" "e3/e4" ok) ("CT" "CT" "e2/e4" ambiguous) ("CT" "TT" "e1/e2" unusual)
                 ("CC" "CC" "e4/e4" ok) ("CC" "CT" "e1/e4" unusual) ("CC" "TT" "e1/e1" unusual))))
    (dolist (c cases)
      (let ((r (genetics-apoe-interpret (nth 0 c) (nth 1 c))))
        (should (equal (nth 2 c) (plist-get r :diplotype)))
        (should (eq (nth 3 c) (plist-get r :status)))
        ;; allele order must not matter
        (should (equal (nth 2 c)
                       (plist-get (genetics-apoe-interpret (reverse (nth 0 c)) (reverse (nth 1 c)))
                                  :diplotype)))))
    (let ((d (plist-get (genetics-apoe-interpret "CT" "CT") :description)))
      (should (string-match-p "e2/e4" d))
      (should (string-match-p "e1/e3" d))
      (should (string-match-p "phase" d)))))

(defun genetics-annotate-test--reverse (s) "Reverse S." (concat (reverse (append s nil))))

(ert-deftest genetics-annotate-test-apoe-edge-cases ()
  (should (eq 'no-call (plist-get (genetics-apoe-interpret "--" "CC") :status)))
  (should (eq 'no-call (plist-get (genetics-apoe-interpret "TT" "--") :status)))
  (should (eq 'strand-flip (plist-get (genetics-apoe-interpret "AG" "GA") :status)))
  (should (null (plist-get (genetics-apoe-interpret "AG" "GA") :diplotype)))
  (should (eq 'incomplete (plist-get (genetics-apoe-interpret "T" "C") :status)))
  (should (eq 'incomplete (plist-get (genetics-apoe-interpret "TA" "CC") :status))))

(ert-deftest genetics-annotate-test-apoe-for-kit-and-annotated-snps ()
  (genetics-test-with-env
    (let ((kit (genetics-test-load "23andme-sample.txt")))
      (should (equal "e2/e4" (plist-get (genetics-apoe-for-kit kit) :diplotype)))
      (let ((hits (genetics-annotated-snps kit)))
        (should (= 8 (length hits)))
        (should (equal "APOE" (genetics-annotation-gene (car (car hits)))))))
    (let ((kit (genetics-test-kit-from-rows '(("rs1" "1" 1 "AA")))))
      (should (null (genetics-apoe-for-kit kit))))))

(provide 'genetics-annotate-test)
;;; genetics-annotate-test.el ends here
