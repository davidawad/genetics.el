;;; genetics-vcf-test.el --- VCF parsing, gzip and offset-index mode -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(ert-deftest genetics-vcf-test-gt-conversion ()
  (should (equal "AA" (genetics--gt-genotype "0/0" "A" '("G"))))
  (should (equal "AG" (genetics--gt-genotype "0/1" "A" '("G"))))
  (should (equal "GA" (genetics--gt-genotype "1|0" "A" '("G"))))
  (should (equal "GG" (genetics--gt-genotype "1|1" "A" '("G"))))
  (should (equal "GT" (genetics--gt-genotype "1/2" "A" '("G" "T"))))
  (should (equal "AT" (genetics--gt-genotype "0/2" "A" '("G" "T"))))
  (should (equal "--" (genetics--gt-genotype "./." "A" '("G"))))
  (should (equal "--" (genetics--gt-genotype "./1" "A" '("G"))))
  (should (equal "G" (genetics--gt-genotype "1" "A" '("G"))))
  (should (equal "AT/A" (genetics--gt-genotype "0/1" "AT" '("A"))))
  (should (equal "ag" (downcase (genetics--gt-genotype "0/1" "a" '("g")))))
  (should-error (genetics--gt-genotype "0/3" "A" '("G")) :type 'genetics-parse-error)
  (should-error (genetics--gt-genotype "x/1" "A" '("G")) :type 'genetics-parse-error))

(ert-deftest genetics-vcf-test-eager-parse ()
  (genetics-test-with-env
    (let ((kit (genetics-parse-file (genetics-test-fixture "sample.vcf"))))
      (should (eq 'vcf (genetics-kit-format kit)))
      (should (equal "37" (genetics-kit-build kit)))
      (should (equal "SAMPLE1" (genetics-kit-sample kit)))
      (should-not (genetics-kit-lazy kit))
      (should (= 13 (genetics-kit-snp-count kit)))
      (cl-flet ((gt (id) (genetics-snp-genotype (genetics-kit-get kit id))))
        (should (equal "AA" (gt "rs4477212")))
        (should (equal "AG" (gt "rs3094315")))        ; phased 0|1
        (should (equal "GG" (gt "rs3131972")))        ; no DP field
        (should (equal "CT" (gt "1:900000")))         ; "." id keyed chrom:pos
        (should (equal "GT" (gt "rs9999")))           ; multi-allelic 1/2
        (should (equal "AT" (gt "rs9998")))           ; 0/2
        (should (equal "AT/A" (gt "rs9997")))         ; indel
        (should (equal "--" (gt "rs9996")))
        (should (equal "--" (gt "rs9995")))
        (should (equal "G" (gt "rs5001"))))           ; haploid
      ;; chr prefix normalized, chrM -> MT
      (should (equal "1" (genetics-snp-chrom (genetics-kit-get kit "rs4477212"))))
      (should (equal "MT" (genetics-snp-chrom (genetics-kit-get kit "rs7001"))))
      (should (equal '("1" "19" "X" "MT") (genetics-kit-chromosomes kit)))
      (should (equal "A" (genetics-snp-ref (genetics-kit-get kit "rs4477212"))))
      (should (equal "G,T" (genetics-snp-alt (genetics-kit-get kit "rs9999"))))
      (should (string-match-p "REF" (genetics-kit-strand-note kit))))))

(ert-deftest genetics-vcf-test-build-38-and-sample-selection ()
  (genetics-test-with-env
    (let ((file (make-temp-file "genetics-test-" nil ".vcf")))
      (unwind-protect
          (progn
            (with-temp-file file
              (insert "##fileformat=VCFv4.2\n##reference=GRCh38\n"
                      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\n"
                      "1\t100\trs1\tA\tG\t.\t.\t.\tGT\t0/0\t1/1\n"))
            (let ((kit (genetics-parse-file file)))
              (should (equal "38" (genetics-kit-build kit)))
              (should (equal "A" (genetics-kit-sample kit)))
              (should (equal "AA" (genetics-snp-genotype (genetics-kit-get kit "rs1")))))
            (let ((kit (genetics-parse-file file :sample "B")))
              (should (equal "B" (genetics-kit-sample kit)))
              (should (equal "GG" (genetics-snp-genotype (genetics-kit-get kit "rs1")))))
            (should-error (genetics-parse-file file :sample "Z") :type 'genetics-parse-error))
        (delete-file file)))))

(ert-deftest genetics-vcf-test-no-sample-column ()
  (genetics-test-with-env
    (let ((file (make-temp-file "genetics-test-" nil ".vcf")))
      (unwind-protect
          (progn
            (with-temp-file file
              (insert "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                      "1\t100\trs1\tA\tG\t.\t.\t.\n"))
            (should-error (genetics-parse-file file) :type 'genetics-parse-error))
        (delete-file file)))))

(ert-deftest genetics-vcf-test-gzip ()
  (genetics-test-with-env
    (let* ((gz (genetics-test-gzip (genetics-test-fixture "sample.vcf")))
           (plain (genetics-parse-file (genetics-test-fixture "sample.vcf"))))
      (unwind-protect
          (let ((kit (genetics-parse-file gz)))
            (should (eq 'vcf (genetics-kit-format kit)))
            (should-not (genetics-kit-lazy kit))
            (should (equal (genetics-test-snps plain) (genetics-test-snps kit)))
            (should (equal "SAMPLE1" (genetics-kit-sample kit))))
        (delete-file gz)))))

(ert-deftest genetics-vcf-test-gzip-not-a-vcf ()
  (genetics-test-with-env
    (let* ((gz (genetics-test-gzip (genetics-test-fixture "23andme-sample.txt"))))
      (unwind-protect
          (should-error (genetics-parse-file gz) :type 'genetics-unknown-format)
        (delete-file gz)))))

(ert-deftest genetics-vcf-test-gzip-missing-binary ()
  (genetics-test-with-env
    (let ((gz (genetics-test-gzip (genetics-test-fixture "sample.vcf"))))
      (unwind-protect
          (cl-letf (((symbol-function 'executable-find) (lambda (&rest _) nil)))
            (should-error (genetics-parse-file gz) :type 'genetics-gzip-error))
        (delete-file gz)))))

(defun genetics-vcf-test--all (kit)
  "Return the list of SNPs of KIT in walk order."
  (let (acc)
    (genetics-kit-map-snps kit (lambda (s) (push s acc) nil))
    (nreverse acc)))

(ert-deftest genetics-vcf-test-offset-index-mode ()
  (genetics-test-with-env
    (let* ((eager (genetics-parse-file (genetics-test-fixture "sample.vcf")))
           (genetics-vcf-eager-limit 10)
           (genetics-chunk-size 64)     ; many chunk boundaries inside lines
           (lazy (genetics-parse-file (genetics-test-fixture "sample.vcf"))))
      (should (genetics-kit-lazy lazy))
      (should-not (genetics-kit-table lazy))
      (should (= 13 (genetics-kit-snp-count lazy)))
      (should (equal "37" (genetics-kit-build lazy)))
      (should (equal (genetics-kit-chromosomes eager) (genetics-kit-chromosomes lazy)))
      (should (equal (plist-get (genetics-kit-stats eager) :chrom-counts)
                     (plist-get (genetics-kit-stats lazy) :chrom-counts)))
      (maphash (lambda (id snp)
                 (should (equal snp (genetics-kit-get lazy id))))
               (genetics-kit-table eager))
      (should (null (genetics-kit-get lazy "rs-nope")))
      (should (equal "CT" (genetics-snp-genotype (genetics-kit-get lazy "1:900000"))))
      ;; chromosome ranges are byte offsets and ordered
      (let ((r (cdr (assoc "19" (genetics-kit-ranges lazy)))))
        (should (< (car r) (cdr r))))
      ;; walking yields the same records in the same order
      (should (equal (genetics-vcf-test--all eager) (genetics-vcf-test--all lazy)))
      ;; per-chromosome walk and early stop
      (let (ids)
        (genetics-kit-map-snps lazy (lambda (s) (push (genetics-snp-rsid s) ids) nil) "19")
        (should (equal '("rs7412" "rs429358") ids)))
      (let ((n 0))
        (genetics-kit-map-snps lazy (lambda (_s) (cl-incf n) (when (= n 3) 'stop)))
        (should (= 3 n))))))

(ert-deftest genetics-vcf-test-offset-index-gz-uses-cache-copy ()
  (genetics-test-with-env
    (let* ((gz (genetics-test-gzip (genetics-test-fixture "sample.vcf")))
           (genetics-vcf-eager-limit 10))
      (unwind-protect
          (let ((kit (genetics-parse-file gz)))
            (should (genetics-kit-lazy kit))
            (should (string-prefix-p (file-truename genetics-cache-directory)
                                     (file-truename (genetics-kit-data-file kit))))
            (should (file-exists-p (genetics-kit-data-file kit)))
            (should (equal "AG" (genetics-snp-genotype (genetics-kit-get kit "rs3094315"))))
            ;; second open reuses the decompressed copy
            (let ((mtime (file-attribute-modification-time
                          (file-attributes (genetics-kit-data-file kit)))))
              (should (equal mtime (file-attribute-modification-time
                                    (file-attributes
                                     (genetics-kit-data-file (genetics-parse-file gz))))))))
        (delete-file gz)))))

(ert-deftest genetics-vcf-test-lazy-kit-browse-and-compare ()
  (genetics-test-with-env
    (let* ((eager (genetics-test-load "sample.vcf"))
           (genetics-vcf-eager-limit 10)
           (lazy (genetics-parse-file (genetics-test-fixture "sample.vcf"))))
      (should (= (length (car (genetics-browse-rows eager '(:chrom "1"))))
                 (length (car (genetics-browse-rows lazy '(:chrom "1"))))))
      (let ((r (genetics-compare-kits eager lazy)))
        (should (= 13 (plist-get r :overlap)))
        (should (= 100.0 (plist-get r :concordance)))))))

(provide 'genetics-vcf-test)
;;; genetics-vcf-test.el ends here
