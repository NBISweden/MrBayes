dist_pdf_DATA=	Manual_MrBayes_v3.2.pdf

CLEANFILES=	\
	Manual_MrBayes_v3.2.aux \
	Manual_MrBayes_v3.2.bbl \
	Manual_MrBayes_v3.2.blg \
	Manual_MrBayes_v3.2.fdb_latexmk \
	Manual_MrBayes_v3.2.fls \
	Manual_MrBayes_v3.2.log \
	Manual_MrBayes_v3.2.toc

DISTCLEANFILES=	Manual_MrBayes_v3.2.pdf

if HAVE_LATEX

envvars=	TEXINPUTS=$(srcdir): BIBINPUTS=$(srcdir): BSTINPUTS=$(srcdir):

if HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf:
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/Manual_MrBayes_v3.2

else # !HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf: \
    Manual_MrBayes_v3.2.tex \
    bayes.bib \
    Appendix_Fig1.pdf \
    Appendix_Fig2.pdf \
    Appendix_Fig3.pdf
	env $(envvars) $(PDFLATEX) Manual_MrBayes_v3.2
	env $(envvars) $(BIBTEX) Manual_MrBayes_v3.2
	env $(envvars) $(PDFLATEX) Manual_MrBayes_v3.2
	env $(envvars) $(PDFLATEX) Manual_MrBayes_v3.2

endif
else # !HAVE_LATEX

Manual_MrBayes_v3.2.pdf:
	echo 'Can not build the manual PDF since LaTeX was not detected' >&2
	exit 0

endif
