all:	Manual_MrBayes_v3.2.pdf

EXTRA_DIST=	Manual_MrBayes_v3.2.tex \
		bayes.bib \
		Appendix_Fig1.pdf \
		Appendix_Fig2.pdf \
		Appendix_Fig3.pdf \
		StatisticalMethods_Chapter.tex \
		sysbio.bst \
		wiki2latex.rb

CLEANFILES=	Manual_MrBayes_v3.2.aux \
		Manual_MrBayes_v3.2.bbl \
		Manual_MrBayes_v3.2.blg \
		Manual_MrBayes_v3.2.log \
		Manual_MrBayes_v3.2.toc

DISTCLEANFILES=	Manual_MrBayes_v3.2.pdf

Manual_MrBayes_v3.2.pdf: \
    Manual_MrBayes_v3.2.tex \
    bayes.bib \
    Appendix_Fig1.pdf \
    Appendix_Fig2.pdf \
    Appendix_Fig3.pdf
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/Manual_MrBayes_v3.2
	BSTINPUTS=$(srcdir) BIBINPUTS=$(srcdir) $(BIBTEX) Manual_MrBayes_v3.2
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/Manual_MrBayes_v3.2
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/Manual_MrBayes_v3.2
