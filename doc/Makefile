# Build the lassip manual.
#
#   make            build ../LASSI-Plus-Manual.pdf
#   make options    regenerate the option reference and the version from the binary
#   make check      fail if either has drifted from the binary, or the README's
#                   --help block has
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

.PHONY: all options check check-options check-version check-readme clean

all: $(PDF)

# Run until the cross-references stop moving rather than a fixed number of
# times. Two passes are not enough here: the second pass inserts the table of
# contents, which pushes every section onto the following page, so the page
# numbers the first pass recorded are all off by one until a third pass.
$(PDF): MANUAL.tex options.tex version.tex Makefile
	$(LATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex
	@n=0; while grep -q "Rerun to get" MANUAL.log && [ $$n -lt 4 ]; do \
		echo "  cross-references moved, rerunning $(LATEX)"; \
		$(LATEX) -interaction=nonstopmode -halt-on-error MANUAL.tex >/dev/null || exit 1; \
		n=$$((n+1)); \
	done; \
	if grep -q "Rerun to get" MANUAL.log; then \
		echo "$(LATEX): cross-references still moving after $$n reruns" >&2; exit 1; \
	fi
	mv MANUAL.pdf $(PDF)

options: options.tex version.tex

options.tex: gen-options.py $(LASSIP)
	$(PYTHON) gen-options.py $(LASSIP) > $@.tmp
	mv $@.tmp $@

# The version on the title page comes from the binary too. It used to be typed
# into MANUAL.tex, which is a second place for it to live and drift from
# src/lassip-cli.h -- the failure this whole arrangement exists to prevent.
version.tex: $(LASSIP)
	@printf '%% Generated from `$(LASSIP) --version`. Do not edit.\n' > $@.tmp
	@printf '\\newcommand{\\lassipversion}{%s}\n' \
		"$$($(LASSIP) --version | awk '{print $$2}' | tr -d v)" >> $@.tmp
	mv $@.tmp $@

check: check-options check-version check-readme

check-options:
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

check-version:
	@$(MAKE) --no-print-directory version.tex >/dev/null
	@v=$$($(LASSIP) --version | awk '{print $$2}' | tr -d v); \
	if grep -q "{$$v}" version.tex; then \
		echo "the manual's version matches the binary ($$v)"; \
	else \
		echo "version.tex does not name $$v; run 'make options'" >&2; exit 1; \
	fi

# The README embeds a copy of --help between the usage line and the change log.
# Nothing regenerates it automatically, so it is the piece most likely to drift.
check-readme:
	@$(LASSIP) --help | awk 'NF{p=NR} {l[NR]=$$0} END{for(i=1;i<=p;i++) print l[i]}' > readme.new
	@n=$$(wc -l < readme.new | tr -d ' '); \
	awk -v n=$$n '/^Usage: lassip --vcf/{f=1} f{print; if(++c==n) exit}' ../README.md > readme.cur
	@if cmp -s readme.cur readme.new; then \
		rm -f readme.cur readme.new; \
		echo "the README's --help block is up to date"; \
	else \
		echo "the README's --help block is out of date; replace it with the output of '$(LASSIP) --help'" >&2; \
		diff -u readme.cur readme.new | head -40 >&2; \
		rm -f readme.cur readme.new; \
		exit 1; \
	fi

clean:
	rm -f MANUAL.aux MANUAL.log MANUAL.out MANUAL.toc MANUAL.pdf options.new options.tex.tmp version.tex.tmp readme.cur readme.new
