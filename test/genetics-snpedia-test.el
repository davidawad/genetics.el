;;; genetics-snpedia-test.el --- SNPedia (mocked) and privacy checks -*- lexical-binding: t; -*-

;;; Code:

(require 'genetics-test-util)
(require 'genetics-snpedia)

(defconst genetics-snpedia-test--json
  "{\"parse\":{\"title\":\"Rs429358\",\"pageid\":1,\"wikitext\":{\"*\":\"{{Rsnum\\n|rsid=429358\\n|Gene=APOE\\n|Orientation=plus\\n|Summary=A fake test summary\\n}}\"}}}"
  "Fake SNPedia API answer.")

(defun genetics-snpedia-test--response (status body)
  "Return a buffer holding an HTTP response with STATUS and BODY."
  (let ((buf (generate-new-buffer " *fake-url*")))
    (with-current-buffer buf
      (set-buffer-multibyte nil)
      (insert (format "HTTP/1.1 %s\r\nContent-Type: application/json\r\n\r\n%s" status body)))
    buf))

(defmacro genetics-snpedia-test-with-mock (status body urls &rest forms)
  "Run FORMS with `url-retrieve-synchronously' mocked; record URLs in URLS."
  (declare (indent 3))
  `(cl-letf (((symbol-function 'url-retrieve-synchronously)
              (lambda (url &rest _)
                (push url ,urls)
                (genetics-snpedia-test--response ,status ,body)))
             ((symbol-function 'url-retrieve)
              (lambda (&rest _) (error "Unexpected async network access"))))
     ,@forms))

(ert-deftest genetics-snpedia-test-disabled-refuses ()
  (genetics-test-with-env
    (let ((urls nil))
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        (should-error (genetics-snpedia-fetch "rs429358") :type 'genetics-snpedia-disabled)
        (should-error (genetics-snpedia-summary "rs429358") :type 'genetics-snpedia-error)
        (should (null urls)))
      (should-not (default-value 'genetics-snpedia-enabled)))))

(ert-deftest genetics-snpedia-test-request-contains-only-rsid ()
  (genetics-test-with-env
    (let ((genetics-snpedia-enabled t) (genetics-snpedia-confirmed t) (urls nil))
      (genetics-test-load "23andme-sample.txt")
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        (let ((s (genetics-snpedia-summary "rs429358")))
          (should (equal "APOE" (plist-get s :gene)))
          (should (equal "plus" (plist-get s :orientation)))
          (should (equal "A fake test summary" (plist-get s :summary))))
        (should (= 1 (length urls)))
        (let ((url (car urls)))
          (should (equal "https://bots.snpedia.com/api.php?action=parse&page=Rs429358&format=json&prop=wikitext" url))
          (should (string-prefix-p "https://bots.snpedia.com/" url))
          ;; nothing about the user's data: no genotype, position, file or kit name
          (let ((case-fold-search nil))
            (dolist (secret '("45411941" "23andme-sample" "fixtures" "heterozygous"))
              (should-not (string-match-p (regexp-quote secret) url)))))))))

(ert-deftest genetics-snpedia-test-disk-cache ()
  (genetics-test-with-env
    (let ((genetics-snpedia-enabled t) (genetics-snpedia-confirmed t) (urls nil))
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        (genetics-snpedia-fetch "rs429358")
        (genetics-snpedia-fetch "rs429358")
        (should (= 1 (length urls)))
        (should (file-exists-p (expand-file-name "snpedia/rs429358.json" genetics-cache-directory)))))))

(ert-deftest genetics-snpedia-test-rejects-non-rsids ()
  (genetics-test-with-env
    (let ((genetics-snpedia-enabled t) (genetics-snpedia-confirmed t) (urls nil))
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        (dolist (bad '("1:900000" "rs12; DROP" "i3000001" "" "rs" "../rs1"))
          (should-error (genetics-snpedia-fetch bad) :type 'genetics-snpedia-error))
        (should (null urls))))))

(ert-deftest genetics-snpedia-test-first-use-confirmation ()
  (genetics-test-with-env
    (let ((genetics-snpedia-enabled t) (urls nil) (asked 0))
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        (cl-letf (((symbol-function 'y-or-n-p) (lambda (&rest _) (cl-incf asked) nil)))
          (should-error (genetics-snpedia-fetch "rs429358") :type 'genetics-snpedia-declined)
          (should (= 1 asked))
          (should (null urls)))
        (cl-letf (((symbol-function 'y-or-n-p) (lambda (&rest _) (cl-incf asked) t)))
          (genetics-snpedia-fetch "rs429358")
          (genetics-snpedia-fetch "rs7412")   ; asks only once per session
          (should (= 2 asked))
          (should (= 2 (length urls))))))))

(ert-deftest genetics-snpedia-test-http-errors ()
  (genetics-test-with-env
    (let ((genetics-snpedia-enabled t) (genetics-snpedia-confirmed t) (urls nil))
      (genetics-snpedia-test-with-mock "404 Not Found" "nope" urls
        (should-error (genetics-snpedia-fetch "rs1") :type 'genetics-snpedia-error))
      (genetics-snpedia-test-with-mock "200 OK" "not json" urls
        (should-error (genetics-snpedia-fetch "rs2") :type 'genetics-snpedia-error))
      (genetics-snpedia-test-with-mock "200 OK" "{\"error\":{\"code\":\"missingtitle\"}}" urls
        (should-error (genetics-snpedia-fetch "rs3") :type 'genetics-snpedia-error))
      (cl-letf (((symbol-function 'url-retrieve-synchronously) (lambda (&rest _) nil)))
        (should-error (genetics-snpedia-fetch "rs4") :type 'genetics-snpedia-error))
      (should-not (file-exists-p (expand-file-name "snpedia/rs1.json" genetics-cache-directory))))))

(ert-deftest genetics-snpedia-test-lookup-integration ()
  (genetics-test-with-env
    (genetics-test-load "23andme-sample.txt")
    (let ((urls nil))
      (genetics-snpedia-test-with-mock "200 OK" genetics-snpedia-test--json urls
        ;; disabled: nothing fetched, nothing shown
        (should-not (string-match-p "^SNPedia$" (genetics-test-buffer-text (genetics-lookup "rs429358"))))
        (should (null urls))
        (let ((genetics-snpedia-enabled t) (genetics-snpedia-confirmed t))
          (let ((text (genetics-test-buffer-text (genetics-lookup "rs429358"))))
            (should (string-match-p "^SNPedia$" text))
            (should (string-match-p "A fake test summary" text))
            (should (= 1 (length urls))))
          ;; a failing lookup is reported in the buffer, not swallowed
          (cl-letf (((symbol-function 'url-retrieve-synchronously) (lambda (&rest _) nil)))
            (should (string-match-p "SNPedia lookup failed"
                                    (genetics-test-buffer-text (genetics-lookup "rs7412"))))))))))

(ert-deftest genetics-privacy-test-network-code-only-in-snpedia ()
  "Network primitives may only appear in genetics-snpedia.el."
  (let* ((root (genetics-test-root))
         (files (directory-files root t "\\`genetics.*\\.el\\'"))
         (rx "url-retrieve\\|make-network-process\\|open-network-stream\\|url-copy-file\\|url-insert-file-contents\\|make-process\\|start-process"))
    (should (> (length files) 8))
    (dolist (f files)
      (let ((text (with-temp-buffer (insert-file-contents f) (buffer-string)))
            (name (file-name-nondirectory f)))
        (if (equal name "genetics-snpedia.el")
            (should (string-match-p "url-retrieve-synchronously" text))
          (should-not (string-match-p rx text)))))
    ;; only snpedia pulls in the url library
    (dolist (f files)
      (unless (equal (file-name-nondirectory f) "genetics-snpedia.el")
        (should-not (string-match-p "(require 'url" (with-temp-buffer (insert-file-contents f) (buffer-string))))))))

(ert-deftest genetics-privacy-test-loading-package-does-not-load-url-snpedia ()
  (should (featurep 'genetics))
  (should-not (default-value 'genetics-snpedia-enabled)))

(provide 'genetics-snpedia-test)
;;; genetics-snpedia-test.el ends here
