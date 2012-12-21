# Makefile for documentation
#
# D. M. Luchtenburg
# Princeton University
#

all: doc # manual.pdf TODO: write manual

.PHONY: doc clean distclean
.SUFFIXES: .tex .pdf
RMFLAGS = -rf

doc:
	@doxygen

.tex.pdf:
	@-(pdflatex $<; pdflatex $<; pdflatex $<)

clean:
	@-$(RM) *.aux *.toc *.blg *.fdb_latexmk *.log *.out *.pdfsync *~
	@-$(RM) -rf ./html ./latex

distclean: clean
	@-$(RM) $(RMFLAGS) html latex
