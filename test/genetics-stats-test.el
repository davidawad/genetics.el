;;; genetics-stats-test.el --- sex inference and contig folding -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)

(defun genetics-stats-test--male-array-rows (build)
  "Synthetic male array rows on BUILD (\"37\" or \"38\").
Non-PAR X is homozygous with 2% heterozygous genotyping errors, PAR1 and
PAR2 are diploid and heterozygous, and 87% of Y sites are called."
  (let* ((par (cdr (assoc build genetics-par-regions)))
         (p1 (car par)) (p2 (cadr par)))
    (append
     ;; 100 non-PAR X sites, 2 heterozygous
     (cl-loop for i from 1 to 100
              collect (list (format "rsx%d" i) "X" (+ 5000000 (* i 1000))
                            (if (<= i 2) "AG" "AA")))
     ;; 40 PAR sites, half heterozygous (diploid in males)
     (cl-loop for i from 1 to 20
              collect (list (format "rsp%d" i) "X" (+ (car p1) (* i 100))
                            (if (cl-evenp i) "CT" "CC")))
     (cl-loop for i from 1 to 20
              collect (list (format "rsq%d" i) "X" (+ (car p2) (* i 100))
                            (if (cl-evenp i) "CT" "TT")))
     ;; 100 Y sites, 87 called
     (cl-loop for i from 1 to 100
              collect (list (format "rsy%d" i) "Y" (* i 1000)
                            (if (<= i 87) "A" "--"))))))

(ert-deftest genetics-stats-test-par-regions ()
  (should (genetics-par-p 60001 "37"))
  (should (genetics-par-p 2699520 "37"))
  (should-not (genetics-par-p 2699521 "37"))
  (should (genetics-par-p 155000000 "37"))
  (should-not (genetics-par-p 2781480 "38"))
  (should-not (genetics-par-p 10001 "37"))
  (should (genetics-par-p 10001 "38"))
  (should (genetics-par-p 156030895 "38"))
  (should-not (genetics-par-p 156030896 "38"))
  ;; unknown build: both builds' regions are excluded
  (should (genetics-par-p 10001 nil))
  (should (genetics-par-p 60001 nil)))

(ert-deftest genetics-stats-test-male-array-with-par-is-male ()
  "A male kit with ~87% Y calls and heterozygous PAR sites is male."
  (genetics-test-with-env
    (let* ((kit (genetics-test-kit-from-rows
                 (genetics-stats-test--male-array-rows "37")))
           (s (genetics-kit-stats kit))
           (e (genetics-sex-evidence kit)))
      (should (= 40 (plist-get s :x-par)))
      (should (= 2 (plist-get s :x-het)))
      (should (< (abs (- 0.02 (plist-get e :x-het-rate))) 1e-9))
      (should (< (abs (- 0.87 (plist-get e :y-called-rate))) 1e-9))
      (should (eq 'male (genetics-infer-sex kit)))
      (should (string-match-p "\\`male (non-PAR X heterozygosity 2.00%, Y called 87%)"
                              (genetics-sex-description kit))))
    ;; the same kit declared GRCh38 uses the GRCh38 PAR coordinates
    (let ((kit (genetics-test-kit-from-rows
                (genetics-stats-test--male-array-rows "38")
                "# We are using reference human assembly build 38 (GRCh38).\n")))
      (should (equal "38" (genetics-kit-build kit)))
      (should (= 40 (plist-get (genetics-kit-stats kit) :x-par)))
      (should (eq 'male (genetics-infer-sex kit))))))

(ert-deftest genetics-stats-test-female-and-conflicts ()
  (genetics-test-with-env
    ;; female: X heterozygous, Y mostly no-call (a few cross-hybridising calls)
    (let ((kit (genetics-test-kit-from-rows
                (append
                 (cl-loop for i from 1 to 100
                          collect (list (format "rsx%d" i) "X" (+ 5000000 i)
                                        (if (zerop (% i 4)) "AG" "GG")))
                 (cl-loop for i from 1 to 100
                          collect (list (format "rsy%d" i) "Y" i
                                        (if (<= i 10) "A" "--")))))))
      (should (eq 'female (genetics-infer-sex kit))))
    ;; heterozygous X and fully called Y (e.g. XXY): uncertain
    (let ((kit (genetics-test-kit-from-rows
                (append
                 (cl-loop for i from 1 to 100
                          collect (list (format "rsx%d" i) "X" (+ 5000000 i)
                                        (if (zerop (% i 4)) "AG" "GG")))
                 (cl-loop for i from 1 to 100
                          collect (list (format "rsy%d" i) "Y" i "A"))))))
      (should (eq 'uncertain (genetics-infer-sex kit))))))

(ert-deftest genetics-stats-test-fold-other-contigs ()
  (should (equal (genetics-fold-chrom-counts
                  '(("1" . 10) ("X" . 3) ("MT" . 1)
                    ("1_KI270706V1_RANDOM" . 4) ("HLA-A*01:01:01:01" . 2)
                    ("UN_GL000220V1" . 1)))
                 '(("1" . 10) ("X" . 3) ("MT" . 1) ("other contigs (3)" . 7))))
  (should (equal (genetics-fold-chrom-counts '(("1" . 5) ("other_contigs" . 9)))
                 '(("1" . 5) ("other contigs" . 9))))
  (should (equal (genetics-fold-chrom-counts '(("1" . 5))) '(("1" . 5)))))

(provide 'genetics-stats-test)
;;; genetics-stats-test.el ends here
