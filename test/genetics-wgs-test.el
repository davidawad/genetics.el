;;; genetics-wgs-test.el --- position-based annotations on rsid-free WGS VCFs -*- lexical-binding: t; -*-

;;; Commentary:

;; wgs-grch38.vcf is synthetic: GRCh38, '.' in every ID, variant sites
;; only, plus alt/decoy/HLA contigs, like a consumer WGS VCF.

;;; Code:

(require 'genetics-test-util)

(defun genetics-wgs-test--load (&rest args)
  "Load the synthetic WGS fixture as variant-only, with extra ARGS."
  (apply #'genetics-test-load "wgs-grch38.vcf" :ref-calls 'absent-means-ref args))

(ert-deftest genetics-wgs-test-curated-coordinates-load ()
  (genetics-test-with-env
    (let ((ann (genetics-annotation "rs429358")))
      (should (equal (genetics-annotation-site ann "38")
                     '(:chrom "19" :pos 44908684 :ref "T")))
      (should (equal (genetics-annotation-site ann "GRCh37")
                     '(:chrom "19" :pos 45411941 :ref "T")))
      (should (cl-some (lambda (u) (string-match-p "ncbi.nlm.nih.gov/snp/rs429358" u))
                       (genetics-annotation-sources ann))))
    ;; every curated entry has both builds, cited
    (maphash (lambda (rsid ann)
               (should (genetics-annotation-site ann "37"))
               (should (genetics-annotation-site ann "38"))
               (should (genetics-annotation-sources ann))
               (should (string-prefix-p "rs" rsid)))
             (genetics-annotations))
    (should (eq (genetics-annotation-at "GRCh38" "1" 11796321)
                (genetics-annotation "rs1801133")))
    (should-not (genetics-annotation-at "GRCh38" "1" 11856378))
    (should (eq (genetics-annotation-at "37" "1" 11856378)
                (genetics-annotation "rs1801133")))))

(ert-deftest genetics-wgs-test-org-coordinates-and-errors ()
  (genetics-test-with-env
    (let* ((anns (genetics-load-annotation-file
                  (expand-file-name "annotations/genetics-example.org"
                                    (genetics-test-root))))
           (mc1r (cl-find "rs1805007" anns :key #'genetics-annotation-rsid
                          :test #'equal)))
      (should (equal (genetics-annotation-site mc1r "38")
                     '(:chrom "16" :pos 89919709 :ref "C"))))
    (let ((file (make-temp-file "genetics-ann" nil ".json")))
      (unwind-protect
          (progn
            (with-temp-file file
              (insert "[{\"rsid\":\"rs1\",\"coordinates\":{\"GRCh99\":{\"chrom\":\"1\",\"pos\":5}}}]"))
            (should-error (genetics-load-annotation-file file)
                          :type 'genetics-annotation-error)
            (with-temp-file file
              (insert "[{\"rsid\":\"rs1\",\"coordinates\":{\"GRCh38\":{\"chrom\":\"1\",\"pos\":-5}}}]"))
            (should-error (genetics-load-annotation-file file)
                          :type 'genetics-annotation-error))
        (delete-file file)))))

(ert-deftest genetics-wgs-test-kit-model ()
  (genetics-test-with-env
    (let ((auto (genetics-parse-file (genetics-test-fixture "wgs-grch38.vcf"))))
      ;; too small for auto-detection: never infer
      (should (eq 'unknown (genetics-kit-ref-calls auto)))
      (should (genetics-kit-rsid-free-p auto)))
    (let ((genetics-wgs-min-records 10))
      (should (eq 'absent-means-ref
                  (genetics-kit-ref-calls
                   (genetics-parse-file (genetics-test-fixture "wgs-grch38.vcf"))))))
    ;; explicit 0/0 calls in the small sample VCF mean it is not variant-only
    (let ((genetics-wgs-min-records 10))
      (should (eq 'unknown (genetics-kit-ref-calls
                            (genetics-parse-file (genetics-test-fixture "sample.vcf"))))))
    (let ((kit (genetics-wgs-test--load)))
      (should (equal "38" (genetics-kit-build kit)))
      (should (eq 'wgs (genetics-kit-assay kit)))
      (should (eq 'male (genetics-infer-sex kit))))))

(ert-deftest genetics-wgs-test-resolve-by-position-and-inferred ()
  (genetics-test-with-env
    (dolist (lazy '(nil t))
      (let* ((kit (genetics-wgs-test--load :lazy lazy))
             (mthfr (genetics-kit-resolve kit "rs1801133"))
             (apoe1 (genetics-kit-resolve kit "rs429358"))
             (apoe2 (genetics-kit-resolve kit "rs7412"))
             (lct (genetics-kit-resolve kit "rs4988235")))
        ;; observed by position
        (should (equal "GA" (genetics-snp-genotype mthfr)))
        (should-not (genetics-snp-inferred-p mthfr))
        (should (equal "CT" (genetics-snp-genotype apoe2)))
        ;; absent from a variant-only VCF: inferred reference, labelled
        (should (equal "TT" (genetics-snp-genotype apoe1)))
        (should (genetics-snp-inferred-p apoe1))
        (should (equal "TT (inferred ref)" (genetics-snp-genotype-label apoe1)))
        (should (equal "GG" (genetics-snp-genotype lct)))
        ;; GRCh38 reference at Factor V Leiden is C
        (should (equal "CC" (genetics-snp-genotype (genetics-kit-resolve kit "rs6025"))))
        ;; a deletion spans the HFE site: never inferred
        (should-not (genetics-kit-resolve kit "rs1800562"))
        ;; chromosome 10 absent altogether: not inferred
        (should-not (genetics-kit-resolve kit "rs4244285"))
        (let ((apoe (genetics-apoe-for-kit kit)))
          (should (equal "e2/e3" (plist-get apoe :diplotype)))
          (should (equal '("rs429358") (plist-get apoe :inferred)))
          (should (string-match-p "inferred homozygous reference" (plist-get apoe :description))))))))

(ert-deftest genetics-wgs-test-no-inference-unless-variant-only ()
  (genetics-test-with-env
    (let ((kit (genetics-test-load "wgs-grch38.vcf")))
      (should (genetics-kit-resolve kit "rs7412"))
      (should-not (genetics-kit-resolve kit "rs429358"))
      (should-not (genetics-apoe-for-kit kit))
      (should (string-match-p "never as homozygous reference"
                              (string-join (genetics-kit-caveats kit) "\n"))))
    ;; arrays list every assayed site: absence is not reference
    (let ((kit (genetics-test-kit-from-rows '(("rs7412" "19" 45412079 "CT")))))
      (should (eq 'explicit (genetics-kit-ref-calls kit)))
      (should-not (genetics-kit-resolve kit "rs429358")))))

(ert-deftest genetics-wgs-test-lazy-block-seek ()
  (genetics-test-with-env
    (let* ((genetics--block-size 2)
           (kit (genetics-wgs-test--load :lazy t)))
      (should (genetics-kit-lazy kit))
      (should (> (length (cdr (assoc "X" (genetics-kit-blocks kit)))) 3))
      (should (equal '(7000000 8000000 9000000)
                     (mapcar #'genetics-snp-pos
                             (genetics-kit-records-in kit "X" 6500000 9000000))))
      (should (equal "AG" (genetics-snp-genotype (genetics-kit-at kit "X" 155800000))))
      (should-not (genetics-kit-at kit "X" 155800001)))))

(ert-deftest genetics-wgs-test-report-summary-lookup ()
  (genetics-test-with-env
    (let* ((kit (genetics-wgs-test--load))
           (report (genetics-report-string kit))
           (summary (genetics-summary-string kit)))
      (dolist (re '("Build: GRCh38" "Assay: wgs (reference calls: absent-means-ref)"
                    "matched by their GRCh38 position"
                    "| rs1801133 | MTHFR | GA | A | 1 |"
                    "| rs429358 | APOE | TT (inferred ref) | C | 0 | Inferred, not observed"
                    "| rs6025 | F5 | CC (inferred ref) | T | 0 |"
                    "Diplotype: e2/e3" "rs429358 inferred homozygous reference"))
        (should (string-match-p (regexp-quote re) report)))
      (should-not (string-match-p "| rs1800562 " report))
      (should (string-match-p "Genotypes marked \"(inferred ref)\" were not observed" report))
      ;; four non-primary contigs folded into one row
      (should (string-match-p "other contigs (4) +4" summary))
      (should-not (string-match-p "HLA\\|EBV\\|KI270706" summary))
      (should (string-match-p "Sex: +male" summary))
      (should (string-match-p "absent = reference, inferred" summary))
      (let ((text (genetics-test-buffer-text (genetics-lookup "rs429358"))))
        (should (string-match-p "TT (inferred ref)  chr19:44908684" text))
        (should (string-match-p "inferred, not observed" text))
        (should (string-match-p "APOE haplotype in wgs-grch38: e2/e3" text)))
      (let ((text (genetics-test-buffer-text (genetics-lookup "rs1801133"))))
        (should (string-match-p "GA  chr1:11796321 .* matched by position" text)))
      ;; browsing: position ids show the gene and count as annotated
      (let ((rows (car (genetics-browse-rows kit '(:annotated t)))))
        (should (equal '("1:11796321" "19:44908822")
                       (mapcar #'genetics-snp-rsid rows))))
      (let ((text (genetics-test-buffer-text (genetics-lookup "19:44908822"))))
        (should (string-match-p "Gene: +APOE" text))))))




;;;; Regressions found in review

(defun genetics-wgs-test--write-vcf (rows &optional format-field)
  "Write a synthetic GRCh38 VCF with ROWS (CHROM POS REF ALT GT); return path.
FORMAT-FIELD is the FORMAT column (default GT, sample value is GT)."
  (let ((file (make-temp-file "genetics-wgs" nil ".vcf")))
    (with-temp-file file
      (insert "##fileformat=VCFv4.2\n##reference=GRCh38\n"
              "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
      (dolist (r rows)
        (insert (format "chr%s\t%d\t.\t%s\t%s\t50\tPASS\t.\t%s\t%s\n"
                        (nth 0 r) (nth 1 r) (nth 2 r) (nth 3 r)
                        (or format-field "GT") (nth 4 r)))))
    file))

(ert-deftest genetics-wgs-test-ref-calls-not-frozen-by-cache ()
  (genetics-test-with-env
    (let ((genetics-use-cache t)
          (file (genetics-test-fixture "wgs-grch38.vcf")))
      (should (eq 'unknown (genetics-kit-ref-calls (genetics-parse-file file))))
      ;; now served from the cache, but the option still decides
      (let ((genetics-vcf-ref-calls 'absent-means-ref))
        (should (eq 'absent-means-ref
                    (genetics-kit-ref-calls (genetics-parse-file file)))))
      (should (eq 'unknown (genetics-kit-ref-calls (genetics-parse-file file)))))))

(ert-deftest genetics-wgs-test-indel-at-site-is-not-the-snp ()
  (genetics-test-with-env
    (let ((file (genetics-wgs-test--write-vcf
                 '((19 44908684 "TG" "T" "0/1")       ; deletion anchored at rs429358
                   (19 44908822 "CA" "C" "0/1")       ; indel, then the SNP, at rs7412
                   (19 44908822 "C" "T" "0/1")))))
      (unwind-protect
          (dolist (lazy '(nil t))
            (let ((kit (genetics-register-kit
                        (genetics-parse-file file :ref-calls 'absent-means-ref :lazy lazy))))
              ;; anchored indel: the SNP is not observed, so it is inferred
              (let ((snp (genetics-kit-resolve kit "rs429358")))
                (should (equal "TT" (genetics-snp-genotype snp)))
                (should (genetics-snp-inferred-p snp)))
              (should (equal "CT" (genetics-snp-genotype (genetics-kit-resolve kit "rs7412"))))
              (should (equal "e2/e3" (plist-get (genetics-apoe-for-kit kit) :diplotype)))
              (unless lazy
                ;; no record is dropped from the eager table
                (should (= 3 (genetics-kit-snp-count kit))))
              (setq genetics-loaded-kits nil)))
        (delete-file file)))))

(ert-deftest genetics-wgs-test-hom-ref-count-uses-sample-gt ()
  (genetics-test-with-env
    (let ((genetics-wgs-min-records 2)
          ;; GT is not the first FORMAT field; 0/0 must still be seen
          (file (genetics-wgs-test--write-vcf
                 '((1 100 "A" "G" "10,0:0/0") (1 200 "C" "T" "5,5:0/1")
                   (1 300 "G" "A" "0,9:1/1"))
                 "AD:GT")))
      (unwind-protect
          (dolist (lazy '(nil t))
            (should (eq 'unknown (genetics-kit-ref-calls
                                  (genetics-parse-file file :lazy lazy)))))
        (delete-file file)))))

(ert-deftest genetics-wgs-test-no-annotation-files ()
  (genetics-test-with-env
    (let ((genetics-annotation-files nil)
          (kit (genetics-wgs-test--load)))
      (should (hash-table-p (genetics-annotations)))
      (should-not (genetics-annotated-snps kit))
      ;; APOE coordinates are built in
      (should (equal "e2/e3" (plist-get (genetics-apoe-for-kit kit) :diplotype))))))

(provide 'genetics-wgs-test)
;;; genetics-wgs-test.el ends here
