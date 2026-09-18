dist_pdf_DATA=	Manual_MrBayes_v3.2.pdf

EXTRA_DIST=	src

DISTCLEANFILES=	Manual_MrBayes_v3.2.pdf

CLEANFILES=	\
	Manual_MrBayes_v3.2.aux \
	Manual_MrBayes_v3.2.bbl \
	Manual_MrBayes_v3.2.blg \
	Manual_MrBayes_v3.2.fdb_latexmk \
	Manual_MrBayes_v3.2.fls \
	Manual_MrBayes_v3.2.log \
	Manual_MrBayes_v3.2.out \
	Manual_MrBayes_v3.2.toc

envvars=	\
	TEXINPUTS=$(srcdir)/src: \
	BIBINPUTS=$(srcdir)/src: \
	BSTINPUTS=$(srcdir)/src:

deps=	\
	$(srcdir)/src/Appendix_Fig1.pdf \
	$(srcdir)/src/Appendix_Fig2.pdf \
	$(srcdir)/src/Appendix_Fig3.pdf \
	$(srcdir)/src/Manual_MrBayes_v3.2.tex \
	$(srcdir)/src/bayes.bib \
	$(srcdir)/src/sysbio.bst

if HAVE_LATEX
if HAVE_LATEXMK

pdf: Manual_MrBayes_v3.2.pdf

Manual_MrBayes_v3.2.pdf: $(deps)
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/src/Manual_MrBayes_v3.2

else # !HAVE_LATEXMK

# The order of pdflatex and bibtex runs below (and the number of runs)
# mirrors what latexmk seems to want to do.
Manual_MrBayes_v3.2.pdf: $(deps)
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Manual_MrBayes_v3.2
	env $(envvars) $(BIBTEX) -terse Manual_MrBayes_v3.2
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Manual_MrBayes_v3.2
	env $(envvars) $(BIBTEX) -terse Manual_MrBayes_v3.2
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Manual_MrBayes_v3.2
	env $(envvars) $(BIBTEX) -terse Manual_MrBayes_v3.2
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Manual_MrBayes_v3.2

endif
else # !HAVE_LATEX

Manual_MrBayes_v3.2.pdf:
	@echo 'Can not build the manual PDF since LaTeX was not detected' >&2
	@exit 0

endif
