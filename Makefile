EMACS ?= emacs
RUN := $(EMACS) -Q --batch -l test/run-tests.el

# All targets go through test/run-tests.el, which also works without make
# (native Windows): emacs -Q --batch -l test/run-tests.el [test|compile|checkdoc|all]
.PHONY: test compile checkdoc lint clean

test:
	$(RUN) test

# Byte-compiles every source with warnings as errors; no .elc is left behind.
compile:
	$(RUN) compile

# Fails (exit 1) when checkdoc reports any warning for a non-test source.
checkdoc:
	$(RUN) checkdoc

lint: compile checkdoc

clean:
	rm -f *.elc
