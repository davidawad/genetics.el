;;; genetics-cache-test.el --- on-disk parse cache -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(ert-deftest genetics-cache-test-round-trip-and-hit ()
  (genetics-test-with-env
    (let ((genetics-use-cache t)
          (calls 0)
          (orig (symbol-function 'genetics--parse-eager))
          (file (genetics-test-fixture "23andme-sample.txt")))
      (cl-letf (((symbol-function 'genetics--parse-eager)
                 (lambda (&rest args) (cl-incf calls) (apply orig args))))
        (let ((first (genetics-parse-file file))
              (second (genetics-parse-file file)))
          (should (= 1 calls))
          (should (directory-files genetics-cache-directory nil "\\.eld\\'"))
          (should (equal (genetics-test-snps first) (genetics-test-snps second)))
          (should (equal (genetics-kit-stats first) (genetics-kit-stats second)))
          (should (equal (genetics-kit-build second) "37"))
          (should (equal (genetics-kit-format second) '23andme))
          (should (equal (genetics-kit-name second) "23andme-sample"))
          (should (equal (genetics-kit-chroms first) (genetics-kit-chroms second))))))))

(ert-deftest genetics-cache-test-invalidated-when-file-changes ()
  (genetics-test-with-env
    (let* ((genetics-use-cache t)
           (file (genetics-test-write-23andme '(("rs1" "1" 10 "AA") ("rs2" "1" 20 "AG")))))
      (unwind-protect
          (progn
            (should (= 2 (genetics-kit-snp-count (genetics-parse-file file))))
            (should (= 2 (genetics-kit-snp-count (genetics-parse-file file))))
            (with-temp-buffer
              (insert "rs3\t1\t30\tGG\n")
              (append-to-file (point-min) (point-max) file))
            (should (= 3 (genetics-kit-snp-count (genetics-parse-file file))))
            ;; same size, different mtime
            (set-file-times file (time-add (current-time) 3600))
            (let ((calls 0) (orig (symbol-function 'genetics--parse-eager)))
              (cl-letf (((symbol-function 'genetics--parse-eager)
                         (lambda (&rest args) (cl-incf calls) (apply orig args))))
                (genetics-parse-file file)
                (should (= 1 calls)))))
        (delete-file file)))))

(ert-deftest genetics-cache-test-corrupt-or-disabled ()
  (genetics-test-with-env
    (let ((file (genetics-test-fixture "23andme-sample.txt")))
      (let ((genetics-use-cache t))
        (genetics-parse-file file)
        (dolist (c (directory-files genetics-cache-directory t "\\.eld\\'"))
          (with-temp-file c (insert "(garbage")))
        (should (= 33 (genetics-kit-snp-count (genetics-parse-file file)))))
      (let ((genetics-use-cache nil)
            (dir (make-temp-file "genetics-test-nocache" t)))
        (unwind-protect
            (let ((genetics-cache-directory dir))
              (genetics-parse-file file)
              (should-not (directory-files dir nil "\\.eld\\'")))
          (delete-directory dir t))))))

(provide 'genetics-cache-test)
;;; genetics-cache-test.el ends here
