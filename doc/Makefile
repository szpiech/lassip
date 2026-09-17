# Build the lassip manual.
#
#   make            build ../LASSI-Plus-Manual.pdf
#   make options    regenerate the option reference from the binary alone
#   make check      fail if options.tex is out of date with the binary
#   make clean      remove the build intermediates
#
# Needs a LaTeX installation (any TeX Live or MacTeX; on macOS pdflatex lives
# in /Library/TeX/texbin, which is not always on PATH) and python3 for the
# option reference. Override either if yours is elsewhere:
#
#   make LATEX=/Library/TeX/texbin/pdflatex
#
# options.tex is generated from `lassip --help`, so the option list in the
# manual describes the binary it was built against. Build the program first;
# the stale v1.1.2 PDF this replaces drifted because nothing tied the two
# together.

LATEX   ?= pdflatex
PYTHON  ?= python3
LASSIP  ?= ../src/lassip
PDF      = ../LASSI-Plus-Manual.pdf

.PHONY: all options check clean

all: $(PDF)

# Two passes: the first writes the .aux and .toc, the second resolves the table
# of contents and the section cross-references against them.
$(PDF): MANUAL.tex options.tex
	$(LATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex
	$(LATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex
	mv MANUAL.pdf $(PDF)

options: options.tex

options.tex: gen-options.py $(LASSIP)
	$(PYTHON) gen-options.py $(LASSIP) > $@.tmp
	mv $@.tmp $@

check:
	@$(PYTHON) gen-options.py $(LASSIP) > options.new
	@if cmp -s options.new options.tex; then \
		rm -f options.new; \
		echo "options.tex is up to date with $(LASSIP) --help"; \
	else \
		echo "options.tex is out of date with $(LASSIP) --help; run 'make options'" >&2; \
		diff -u options.tex options.new | head -40 >&2; \
		rm -f options.new; \
		exit 1; \
	fi

clean:
	rm -f MANUAL.aux MANUAL.log MANUAL.out MANUAL.toc MANUAL.pdf options.new options.tex.tmp
