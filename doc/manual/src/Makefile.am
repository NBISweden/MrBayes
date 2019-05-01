dist_pdf_DATA=	Manual_MrBayes_v3.2.pdf

sources=	\
	$(srcdir)/Appendix_Fig1.pdf \
	$(srcdir)/Appendix_Fig2.pdf \
	$(srcdir)/Appendix_Fig3.pdf \
	$(srcdir)/Manual_MrBayes_v3.2.tex \
	$(srcdir)/StatisticalMethods_Chapter.tex \
	$(srcdir)/bayes.bib \
	$(srcdir)/sysbio.bst

CLEANFILES=	Manual_MrBayes_v3.2.aux \
		Manual_MrBayes_v3.2.bbl \
		Manual_MrBayes_v3.2.blg \
		Manual_MrBayes_v3.2.log \
		Manual_MrBayes_v3.2.toc \
		Manual_MrBayes_v3.2.fdb_latexmk \
		Manual_MrBayes_v3.2.fls \
		$(sources)

DISTCLEANFILES=	Manual_MrBayes_v3.2.pdf

if HAVE_LATEX
if HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf:
	cp $(sources) .
	$(LATEXMK) -quiet -pdf -cd Manual_MrBayes_v3.2

else # !HAVE_LATEXMK

Manual_MrBayes_v3.2.pdf: \
    Manual_MrBayes_v3.2.tex \
    bayes.bib \
    Appendix_Fig1.pdf \
    Appendix_Fig2.pdf \
    Appendix_Fig3.pdf
	cp $(sources) .
	$(PDFLATEX) Manual_MrBayes_v3.2
	$(BIBTEX) Manual_MrBayes_v3.2
	$(PDFLATEX) Manual_MrBayes_v3.2
	$(PDFLATEX) Manual_MrBayes_v3.2

endif
else # !HAVE_LATEX

Manual_MrBayes_v3.2.pdf:
	echo 'Can not build the manual PDF since LaTeX was not detected' >&2
	exit 0

endif
