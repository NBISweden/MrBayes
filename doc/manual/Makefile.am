dist_pdf_DATA=	Manual_MrBayes_v3.2.pdf

EXTRA_DIST=	src

CLEANFILES=	Manual_MrBayes_v3.2.aux \
		Manual_MrBayes_v3.2.bbl \
		Manual_MrBayes_v3.2.blg \
		Manual_MrBayes_v3.2.log \
		Manual_MrBayes_v3.2.toc

DISTCLEANFILES=	Manual_MrBayes_v3.2.pdf

if HAVE_LATEX
if HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf:
	$(LATEXMK) -quiet -pdf -cd $(srcdir)/src/Manual_MrBayes_v3.2
	mv $(srcdir)/src/Manual_MrBayes_v3.2.pdf Manual_MrBayes_v3.2.pdf

else # !HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf: \
    src/Manual_MrBayes_v3.2.tex \
    src/bayes.bib \
    src/Appendix_Fig1.pdf \
    src/Appendix_Fig2.pdf \
    src/Appendix_Fig3.pdf
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/Manual_MrBayes_v3.2
	BSTINPUTS=$(srcdir)/src BIBINPUTS=$(srcdir)/src $(BIBTEX) Manual_MrBayes_v3.2
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/Manual_MrBayes_v3.2
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/Manual_MrBayes_v3.2

endif
else # !HAVE_LATEX

Manual_MrBayes_v3.2.pdf:
	echo 'Can not build the manual PDF since LaTeX was not detected' >&2
	exit 0

endif
