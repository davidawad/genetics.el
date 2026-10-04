;;; genetics-source-test.el --- genome-cli source layer, against a fake genome -*- lexical-binding: t; -*-

;;; Commentary:

;; test/bin/genome replays recorded, synthetic genome/v1 JSON from
;; test/fixtures/genome and logs every argv it receives.

;;; Code:

(require 'genetics-test-util)

(defun genetics-source-test--argv (cmd &rest args)
  "Return the logged form (no executable) of genome argv CMD ARGS."
  (mapconcat #'identity (cdr (apply #'genetics-genome-argv cmd args)) " "))

(ert-deftest genetics-source-test-argv-and-explain ()
  (let ((genetics-genome-executable "genome")
        (genetics-genome-page-size 100))
    (should (equal (genetics-genome-argv 'import "/data/kit.vcf.gz")
                   '("genome" "import" "/data/kit.vcf.gz" "--format" "json")))
    (should (equal (genetics-genome-argv 'summary "k1")
                   '("genome" "summary" "k1" "--format" "json")))
    (should (equal (genetics-genome-argv 'lookup "k1" "rs7412")
                   '("genome" "lookup" "k1" "--rsid" "rs7412" "--format" "json")))
    (should (equal (genetics-genome-argv 'query "k1" "19" 200)
                   '("genome" "export" "k1" "--region" "19" "--limit" "100"
                     "--offset" "200" "--format" "json")))
    (should (equal (genetics-genome-argv 'query "k1" "19" 0 5 9)
                   '("genome" "export" "k1" "--region" "19:5-9"
                     "--limit" "100" "--offset" "0" "--format" "json")))
    (should (equal (genetics-genome-argv 'compare "k1" "k2")
                   '("genome" "compare" "k1" "k2" "--format" "json")))
    (should-error (genetics-genome-argv 'frobnicate) :type 'genetics-genome-error)
    ;; the explain twin is pure and exact
    (should (equal (genetics-source-genome-cli-explain "/data/my kit.txt")
                   "genome import /data/my\\ kit.txt --format json"))))

(ert-deftest genetics-source-test-auto-picks-source ()
  (genetics-test-with-env
    (let ((genetics-source-function #'genetics-source-auto)
          (genetics-genome-executable "genome-cli-surely-not-installed"))
      (should-not (genetics-genome-available-p))
      (let ((kit (genetics-source-open (genetics-test-fixture "23andme-sample.txt"))))
        (should-not (genetics-kit-backend kit))
        (should (= 33 (genetics-kit-snp-count kit)))))
    (genetics-test-with-fake-genome log
      (let ((genetics-source-function #'genetics-source-auto))
        (should (genetics-genome-available-p))
        (should (genetics-kit-backend
                 (genetics-source-open (genetics-test-fixture "wgs-grch38.vcf"))))
        (should (equal "kits --format json" (car (funcall log))))
        (should (string-prefix-p "import " (cadr (funcall log))))))))

(ert-deftest genetics-source-test-open-maps-kit-and-summary ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((file (genetics-test-fixture "wgs-grch38.vcf"))
             (kit (genetics-open file))
             (s (genetics-kit-stats kit)))
        (should (equal (funcall log)
                       (list (concat "import " file " --format json")
                             "summary k1 --format json")))
        (should (equal "k1" (genetics-kit-backend-id kit)))
        (should (eq 'vcf (genetics-kit-format kit)))
        (should (equal "38" (genetics-kit-build kit)))
        (should (eq 'wgs (genetics-kit-assay kit)))
        (should (eq 'absent-means-ref (genetics-kit-ref-calls kit)))
        (should (eq 'none (genetics-kit-has-rsids kit)))
        (should (= 27 (genetics-kit-snp-count kit)))
        (should (= 15 (plist-get s :hom)))
        (should (equal '("1" "2" "6" "19" "X" "Y" "MT")
                       (genetics-kit-chromosomes kit)))
        (should (eq 'male (genetics-infer-sex kit)))
        (let ((text (genetics-test-buffer-text "*genetics: wgs-grch38*")))
          (dolist (re '("Source: +genome-cli (kit k1)" "Build: +GRCh38"
                        "held by genome-cli" "Sex: +male (non-PAR X heterozygosity 9.00%, Y called 100%)"
                        "other contigs +4" "Parsed by genome-cli"
                        "synthetic fixture" "Variant-only VCF: absent sites"
                        "no rsids"))
            (should (string-match-p re text))))))))

(ert-deftest genetics-source-test-lookup-report-apoe ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((kit (genetics-open (genetics-test-fixture "wgs-grch38.vcf")))
             (inferred (genetics-kit-get kit "rs429358")))
        (should (genetics-snp-inferred-p inferred))
        (should (equal "TT" (genetics-snp-genotype inferred)))
        (should (equal "CT" (genetics-snp-genotype (genetics-kit-get kit "rs7412"))))
        ;; call_source "missing" and empty answers are absent
        (should-not (genetics-kit-get kit "rs4244285"))
        (should-not (genetics-kit-get kit "rs6025"))
        ;; memoized: one lookup per rsid
        (genetics-kit-get kit "rs7412")
        (should (= 1 (cl-count "lookup k1 --rsid rs7412 --format json"
                               (funcall log) :test #'equal)))
        (let ((apoe (genetics-apoe-for-kit kit)))
          (should (equal "e2/e3" (plist-get apoe :diplotype)))
          (should (equal '("rs429358") (plist-get apoe :inferred))))
        (let ((report (genetics-report-string kit)))
          (should (string-match-p (regexp-quote "| rs429358 | APOE | TT (inferred ref) |") report))
          (should (string-match-p (regexp-quote "| rs1801133 | MTHFR | GA | A | 1 |") report))
          (should-not (string-match-p "| rs4244285 " report)))
        (should (string-match-p "TT (inferred ref)"
                                (genetics-test-buffer-text (genetics-lookup "rs429358"))))))))

(ert-deftest genetics-source-test-browse-export-paging ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((genetics-genome-page-size 2)
             (kit (genetics-open (genetics-test-fixture "wgs-grch38.vcf")))
             (rows (car (genetics-browse-rows kit '(:chrom "19")))))
        (should (equal '(44908822 44910000 45000000) (mapcar #'genetics-snp-pos rows)))
        (should (equal '("19:44908822" "19:44910000" "19:45000000")
                       (mapcar #'genetics-snp-rsid rows)))
        (should (equal "G,T" (genetics-snp-alt (nth 2 rows))))
        (should (member (genetics-source-test--argv 'query "k1" "19" 2) (funcall log)))
        ;; position ids are annotated by GRCh38 position
        (should (equal '("19:44908822")
                       (mapcar #'genetics-snp-rsid
                               (car (genetics-browse-rows kit '(:chrom "19" :annotated t))))))
        ;; a display limit stops paging early
        (let ((before (length (funcall log)))
              (res (genetics-browse-rows kit '(:chrom "19") 1)))
          (should (= 1 (length (car res))))
          (should (cdr res))
          (should (= 1 (- (length (funcall log)) before))))
        (should (equal "CT" (genetics-snp-genotype (genetics-kit-at kit "19" 44908822))))
        (let ((file (make-temp-file "genetics-export" nil ".csv")))
          (unwind-protect
              (progn
                (should (= 3 (genetics-export-csv file kit '(:chrom "19"))))
                (should (string-match-p "\"19:44908822\",\"19\",\"44908822\",\"CT\",\"heterozygous\",\"APOE\""
                                        (with-temp-buffer (insert-file-contents file) (buffer-string)))))
            (delete-file file)))))))

(ert-deftest genetics-source-test-compare ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome log
      (let* ((a (genetics-open (genetics-test-fixture "wgs-grch38.vcf")))
             (second (make-temp-file "second-" nil ".txt")))
        (unwind-protect
            (let* ((b (genetics-open second))
                   (r (genetics-compare-kits a b)))
              (should (equal "k2" (genetics-kit-backend-id b)))
              (should (eq 'array (genetics-kit-assay b)))
              (should (member "compare k1 k2 --format json" (funcall log)))
              (should (= 10 (plist-get r :overlap)))
              (should (= 9 (plist-get r :concordant)))
              (should (< (abs (- 90.0 (plist-get r :concordance))) 1e-9))
              (should (equal '(("rs7412" "CT" "TT")) (plist-get r :discordant)))
              (let ((text (genetics-test-buffer-text (genetics-compare a b))))
                (should (string-match-p "Concordance: +90.00%" text))
                (should (string-match-p "rs7412 +CT +TT" text))
                (should (string-match-p "lifted from GRCh37" text))))
          (delete-file second))
        ;; mixing sources is refused, not silently wrong
        (let ((native (let ((genetics-source-function #'genetics-source-native))
                        (genetics-source-open (genetics-test-fixture "23andme-sample.txt")))))
          (should-error (genetics-compare-kits a native) :type 'genetics-error))))))

(ert-deftest genetics-source-test-errors ()
  (genetics-test-with-env
    (genetics-test-with-fake-genome _log
      (let ((err (should-error (genetics-source-open "/tmp/missing-kit.txt")
                               :type 'genetics-genome-error)))
        (should (string-match-p "not_found.*no such file" (error-message-string err))))
      (let ((err (should-error (genetics-source-open "/tmp/garbage-kit.txt")
                               :type 'genetics-genome-error)))
        (should (string-match-p "exit 101.*panic: not json" (error-message-string err)))))
    (let ((genetics-genome-executable "genome-cli-surely-not-installed"))
      (should-error (genetics-source-genome-cli (genetics-test-fixture "sample.vcf"))
                    :type 'genetics-genome-missing))))

(ert-deftest genetics-source-test-record-mapping ()
  (let ((snp (genetics-genome-record->snp
              '((rsid) (chrom . "chr19") (pos . 5) (ref . "C") (alt "T")
                (genotype) (zygosity . "no_call") (call_source . "observed")))))
    (should (equal "19:5" (genetics-snp-rsid snp)))
    (should (equal "19" (genetics-snp-chrom snp)))
    (should (equal "--" (genetics-snp-genotype snp)))
    (should-not (genetics-snp-inferred-p snp)))
  (should-not (genetics-genome-record->snp '((call_source . "missing") (chrom . "1") (pos . 1))))
  (should (equal "A" (genetics-snp-genotype
                      (genetics-genome-record->snp
                       '((rsid . "rs1") (chrom . "X") (pos . 9) (genotype . "a")
                         (zygosity . "hemi") (call_source . "observed")))))))


(ert-deftest genetics-source-test-enum-spelling ()
  (should (eq 'absent-means-ref (genetics--genome-symbol "absent_means_ref")))
  (should (eq 'absent-means-ref (genetics--genome-symbol "absent-means-ref")))
  (should (equal '("19:5" "CT" "TT")
                 (genetics--genome-pair '((a . ((chrom . "19") (pos . 5) (genotype . "CT")))
                                          (b . ((chrom . "19") (pos . 5) (genotype . "TT")))))))
  (should (equal '("?" "CT" "TT") (genetics--genome-pair '((a . "CT") (b . "TT"))))))

(provide 'genetics-source-test)
;;; genetics-source-test.el ends here
