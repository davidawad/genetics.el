;;; genetics-parse-test.el --- detection and array-format parsers -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(ert-deftest genetics-parse-test-detect-formats ()
  (should (eq '23andme (genetics-detect-format (genetics-test-fixture "23andme-sample.txt"))))
  (should (eq 'ancestry (genetics-detect-format (genetics-test-fixture "ancestry-sample.txt"))))
  (should (eq 'myheritage (genetics-detect-format (genetics-test-fixture "myheritage-sample.csv"))))
  (should (eq 'ftdna (genetics-detect-format (genetics-test-fixture "ftdna-sample.csv"))))
  (should (eq 'vcf (genetics-detect-format (genetics-test-fixture "sample.vcf")))))

(ert-deftest genetics-parse-test-23andme ()
  (genetics-test-with-env
    (let* ((kit (genetics-parse-file (genetics-test-fixture "23andme-sample.txt")))
           (s (genetics-kit-stats kit)))
      (should (eq '23andme (genetics-kit-format kit)))
      (should (equal "37" (genetics-kit-build kit)))
      (should (equal "unknown" (genetics-kit-chip kit)))
      (should (equal "23andme-sample" (genetics-kit-name kit)))
      (should (= 33 (plist-get s :total)))
      (should (= 6 (plist-get s :nocalls)))
      (should (= 10 (plist-get s :het)))
      (should (= 15 (plist-get s :hom)))
      (should (= 2 (plist-get s :hemi)))
      (should (equal '("1" "2" "3" "6" "10" "19" "X" "Y" "MT")
                     (genetics-kit-chromosomes kit)))
      (should (= 10 (cdr (assoc "X" (plist-get s :chrom-counts)))))
      (let ((snp (genetics-kit-get kit "rs3094315")))
        (should (equal "AG" (genetics-snp-genotype snp)))
        (should (equal "1" (genetics-snp-chrom snp)))
        (should (= 752566 (genetics-snp-pos snp))))
      ;; internal ids and no-calls
      (should (equal "--" (genetics-snp-genotype (genetics-kit-get kit "i3000001"))))
      (should (equal "A" (genetics-snp-genotype (genetics-kit-get kit "rs7001"))))
      (should (string-match-p "forward" (genetics-kit-strand-note kit))))))

(ert-deftest genetics-parse-test-per-chromosome-sorted-by-position ()
  (genetics-test-with-env
    (let ((kit (genetics-parse-file (genetics-test-fixture "23andme-sample.txt"))))
      (let ((positions (mapcar #'genetics-snp-pos
                               (cdr (assoc "1" (mapcar (lambda (c) (cons (car c) (append (cdr c) nil)))
                                                       (genetics-kit-chroms kit)))))))
        (should (equal positions (sort (copy-sequence positions) #'<)))
        (should (= 9 (length positions)))))))

(ert-deftest genetics-parse-test-ancestry ()
  (genetics-test-with-env
    (let ((kit (genetics-parse-file (genetics-test-fixture "ancestry-sample.txt"))))
      (should (eq 'ancestry (genetics-kit-format kit)))
      (should (equal "37" (genetics-kit-build kit)))
      (should (equal "AG" (genetics-snp-genotype (genetics-kit-get kit "rs3094315"))))
      (should (equal "--" (genetics-snp-genotype (genetics-kit-get kit "rs3131972"))))
      (should (equal "X" (genetics-snp-chrom (genetics-kit-get kit "rs5001"))))
      (should (equal "Y" (genetics-snp-chrom (genetics-kit-get kit "rs6001"))))
      (should (equal "XY" (genetics-snp-chrom (genetics-kit-get kit "rs8001"))))
      (should (equal "MT" (genetics-snp-chrom (genetics-kit-get kit "rs7001"))))
      (should (= 12 (genetics-kit-snp-count kit)))
      (should (= 2 (plist-get (genetics-kit-stats kit) :nocalls))))))

(ert-deftest genetics-parse-test-csv-formats ()
  (genetics-test-with-env
    (dolist (case '(("myheritage-sample.csv" . myheritage) ("ftdna-sample.csv" . ftdna)))
      (let ((kit (genetics-parse-file (genetics-test-fixture (car case)))))
        (should (eq (cdr case) (genetics-kit-format kit)))
        (should (equal "37" (genetics-kit-build kit)))
        (should (= 8 (genetics-kit-snp-count kit)))
        (should (equal "AG" (genetics-snp-genotype (genetics-kit-get kit "rs1801133"))))
        (should (equal "--" (genetics-snp-genotype (genetics-kit-get kit "rs3131972"))))
        (should (equal 82154 (genetics-snp-pos (genetics-kit-get kit "rs4477212"))))))
    ;; FTDNA has no comments, so the build is assumed and the note says so
    (should (string-match-p "assumed"
                            (genetics-kit-strand-note
                             (genetics-parse-file (genetics-test-fixture "ftdna-sample.csv")))))))

(ert-deftest genetics-parse-test-build-detection ()
  (should (equal "37" (genetics--detect-build "# reference human assembly build 37 (also known as")))
  (should (equal "37" (genetics--detect-build "#... build 37.1 coordinates")))
  (should (equal "38" (genetics-test--b38)))
  (should (equal "38" (genetics--detect-build "##reference=GRCh38")))
  (should (equal "37" (genetics--detect-build "##reference=hg19")))
  (should (equal "38" (genetics--detect-build "##contig=<ID=1,length=248956422>")))
  (should (null (genetics--detect-build "# nothing here"))))

(defun genetics-test--b38 ()
  "Detect a build 38 comment from a 23andMe file."
  (genetics-test-with-env
    (let ((kit (genetics-test-kit-from-rows
                '(("rs1" "1" 10 "AA"))
                "# We are using reference human assembly build 38.\n")))
      (genetics-kit-build kit))))

(ert-deftest genetics-parse-test-chip-detection ()
  (should (equal "v5" (genetics--detect-chip "" 630000 nil)))
  (should (equal "v4" (genetics--detect-chip "" 570000 nil)))
  (should (equal "v3" (genetics--detect-chip "" 960000 nil)))
  (should (equal "unknown" (genetics--detect-chip "" 1000 nil)))
  (should (equal "v4" (genetics--detect-chip "# chip: v4" 5 nil)))
  (should (equal "v3" (genetics--detect-chip "# chip: v4" 5 "v3")))
  (genetics-test-with-env
    (let ((kit (genetics-test-kit-from-rows '(("rs1" "1" 10 "AA"))
                                            "# 23andMe chip version v5\n")))
      (should (equal "v5" (genetics-kit-chip kit))))
    (let ((file (genetics-test-write-23andme '(("rs1" "1" 10 "AA")))))
      (unwind-protect
          (should (equal "v3" (genetics-kit-chip (genetics-parse-file file :chip "v3"))))
        (delete-file file)))))

(ert-deftest genetics-parse-test-sex-inference ()
  (genetics-test-with-env
    (should (eq 'female (genetics-infer-sex (genetics-parse-file (genetics-test-fixture "23andme-sample.txt")))))
    (let ((male (genetics-test-kit-from-rows
                 (append (cl-loop for i from 1 to 20 collect (list (format "rsx%d" i) "X" (* i 100) (if (cl-evenp i) "A" "G")))
                         (cl-loop for i from 1 to 10 collect (list (format "rsy%d" i) "Y" (* i 100) "A"))))))
      (should (eq 'male (genetics-infer-sex male)))
      (should (string-match-p "^male" (genetics-sex-description male))))
    ;; male with diploid-looking homozygous X calls and Y calls
    (let ((male2 (genetics-test-kit-from-rows
                  (append (cl-loop for i from 1 to 20 collect (list (format "rsx%d" i) "X" (* i 100) "AA"))
                          (cl-loop for i from 1 to 10 collect (list (format "rsy%d" i) "Y" (* i 100) "GG"))))))
      (should (eq 'male (genetics-infer-sex male2))))
    (let ((female (genetics-test-kit-from-rows
                   (append (cl-loop for i from 1 to 20 collect (list (format "rsx%d" i) "X" (* i 100) (if (cl-evenp i) "AG" "CC")))
                           (cl-loop for i from 1 to 10 collect (list (format "rsy%d" i) "Y" (* i 100) "--"))))))
      (should (eq 'female (genetics-infer-sex female))))
    ;; no X data at all: uncertain
    (should (eq 'uncertain (genetics-infer-sex (genetics-test-kit-from-rows '(("rs1" "1" 10 "AA"))))))
    ;; Y called but X heterozygous: uncertain
    (should (eq 'uncertain
                (genetics-infer-sex
                 (genetics-test-kit-from-rows
                  (append (cl-loop for i from 1 to 10 collect (list (format "rsx%d" i) "X" i (if (cl-evenp i) "AG" "CC")))
                          (cl-loop for i from 1 to 10 collect (list (format "rsy%d" i) "Y" i "A")))))))))

(ert-deftest genetics-parse-test-genotype-helpers ()
  (should (eq 'no-call (genetics-zygosity "--")))
  (should (eq 'heterozygous (genetics-zygosity "AG")))
  (should (eq 'homozygous (genetics-zygosity "CC")))
  (should (eq 'hemizygous (genetics-zygosity "A")))
  (should (eq 'heterozygous (genetics-zygosity "DI")))
  (should (equal (genetics-genotype-key "GA") (genetics-genotype-key "AG")))
  (should (equal (genetics-genotype-key "A") (genetics-genotype-key "AA")))
  (should (null (genetics-genotype-key "--")))
  (should (equal "T" (genetics-complement "A")))
  (should (null (genetics-complement "D")))
  (should (equal "MT" (genetics-normalize-chrom "chrM")))
  (should (equal "1" (genetics-normalize-chrom "chr1")))
  (should (equal "X" (genetics-normalize-chrom "23" t)))
  (should (equal "23" (genetics-normalize-chrom "23"))))

(ert-deftest genetics-parse-test-crlf-lines ()
  (genetics-test-with-env
    (let ((file (make-temp-file "genetics-test-" nil ".txt")))
      (unwind-protect
          (progn
            (let ((coding-system-for-write 'binary))
              (with-temp-file file
                (insert "# 23andMe\r\n# rsid\tchromosome\tposition\tgenotype\r\nrs1\t1\t5\tAG\r\nrs2\t1\t6\t--\r\n")))
            (let ((kit (genetics-parse-file file)))
              (should (equal "AG" (genetics-snp-genotype (genetics-kit-get kit "rs1"))))
              (should (equal "--" (genetics-snp-genotype (genetics-kit-get kit "rs2"))))))
        (delete-file file)))))

(defun genetics-parse-test--crlf-copy (name)
  "Return a temp copy of fixture NAME with CRLF line endings."
  (let ((file (make-temp-file "genetics-test-crlf-" nil
                              (concat "." (file-name-extension name)))))
    (with-temp-buffer
      (set-buffer-multibyte nil)
      (insert-file-contents-literally (genetics-test-fixture name))
      (goto-char (point-min))
      (while (search-forward "\n" nil t) (replace-match "\r\n" t t))
      (let ((coding-system-for-write 'no-conversion))
        (write-region nil nil file nil 'silent)))
    file))

(ert-deftest genetics-parse-test-crlf-fixtures ()
  "Every format reads the same from CRLF (Windows) files, eager and lazy."
  (genetics-test-with-env
    (dolist (name '("23andme-sample.txt" "ancestry-sample.txt"
                    "myheritage-sample.csv" "ftdna-sample.csv" "sample.vcf"))
      (let ((file (genetics-parse-test--crlf-copy name))
            (plain (genetics-parse-file (genetics-test-fixture name))))
        (unwind-protect
            (let ((kit (genetics-parse-file file)))
              (should (eq (genetics-kit-format plain) (genetics-kit-format kit)))
              (should (equal (genetics-kit-build plain) (genetics-kit-build kit)))
              (should (equal (genetics-kit-sample plain) (genetics-kit-sample kit)))
              (should (equal (genetics-test-snps plain) (genetics-test-snps kit)))
              (when (eq 'vcf (genetics-kit-format plain))
                (let* ((genetics-vcf-eager-limit 10)
                       (lazy (genetics-parse-file file)))
                  (should (genetics-kit-lazy lazy))
                  (should (equal (genetics-kit-sample plain) (genetics-kit-sample lazy)))
                  (should (equal (genetics-kit-build plain) (genetics-kit-build lazy)))
                  (maphash (lambda (id snp)
                             (should (equal snp (genetics-kit-get lazy id))))
                           (genetics-kit-table plain)))))
          (delete-file file))))))

(ert-deftest genetics-parse-test-errors ()
  (genetics-test-with-env
    (let ((file (make-temp-file "genetics-test-" nil ".txt")))
      (unwind-protect
          (progn
            (with-temp-file file (insert "hello world\n"))
            (should-error (genetics-parse-file file) :type 'genetics-unknown-format)
            (with-temp-file file)
            (should-error (genetics-parse-file file) :type 'genetics-unknown-format)
            (with-temp-file file
              (insert "# 23andMe\nrs1\t1\t5\tAG\nrs2\t1\tabc\tAG\n"))
            (let ((err (should-error (genetics-parse-file file) :type 'genetics-parse-error)))
              (should (string-match-p "line 3" (cadr err))))
            (with-temp-file file (insert "# 23andMe\n# rsid\tchromosome\tposition\tgenotype\n"))
            (should-error (genetics-parse-file file) :type 'genetics-unknown-format))
        (delete-file file)))
    (should-error (genetics-parse-file "/nonexistent/genetics-nope.txt")
                  :type 'genetics-file-error)
    (should (condition-case nil (genetics-parse-file "/nonexistent/x") (genetics-error t)))))

(ert-deftest genetics-parse-test-large-synthetic-file ()
  "A 50k-line file parses within a generous time bound."
  (genetics-test-with-env
    (let ((file (make-temp-file "genetics-test-big-" nil ".txt"))
          (chroms ["1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "X"]))
      (unwind-protect
          (progn
            (with-temp-file file
              (insert "# We are using reference human assembly build 37\n# rsid\tchromosome\tposition\tgenotype\n")
              (dotimes (i 50000)
                (insert (format "rs%d\t%s\t%d\t%s\n" (+ 1000 i) (aref chroms (% i 13))
                                (* 37 (1+ i)) (aref ["AA" "AG" "GG" "--" "CT"] (% i 5))))))
            (let* ((start (float-time))
                   (kit (genetics-parse-file file))
                   (elapsed (- (float-time) start)))
              (should (= 50000 (genetics-kit-snp-count kit)))
              (should (= 10000 (plist-get (genetics-kit-stats kit) :nocalls)))
              (should (equal "CT" (genetics-snp-genotype (genetics-kit-get kit "rs1004"))))
              (should (< elapsed 60.0))
              (message "parsed 50000 lines in %.2fs" elapsed)))
        (delete-file file)))))

(provide 'genetics-parse-test)
;;; genetics-parse-test.el ends here
