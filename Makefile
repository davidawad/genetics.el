EMACS ?= emacs
TESTS := $(wildcard test/*-test.el)
SOURCES := $(filter-out %-test.el,$(wildcard *.el))

.PHONY: test compile checkdoc lint clean

test:
	$(EMACS) -Q --batch -L . -L test $(foreach t,$(TESTS),-l $(t)) -f ert-run-tests-batch-and-exit

compile:
	$(EMACS) -Q --batch -L . --eval '(setq byte-compile-error-on-warn t)' -f batch-byte-compile *.el
	@rm -f *.elc

# Fails (exit 1) when checkdoc reports any warning for a non-test source.
checkdoc:
	$(EMACS) -Q --batch -L . --eval "(progn (require 'checkdoc) (setq checkdoc-diagnostic-buffer \"*genetics-checkdoc*\") (dolist (f command-line-args-left) (with-current-buffer (find-file-noselect f) (checkdoc-current-buffer t))) (let ((b (get-buffer checkdoc-diagnostic-buffer))) (when (and b (with-current-buffer b (goto-char (point-min)) (re-search-forward \"^[^*\\n\\f]\" nil t))) (princ (with-current-buffer b (buffer-string))) (kill-emacs 1))))" $(SOURCES)

lint: compile checkdoc

clean:
	rm -f *.elc
