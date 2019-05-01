dist_pdf_DATA=	MBDating.pdf

CLEANFILES=	MBDating.aux \
		MBDating.bbl \
		MBDating.blg \
		MBDating.log \
		MBDating.out

DISTCLEANFILES=	MBDating.pdf

if HAVE_LATEX

envvars=	TEXINPUTS=$(srcdir): BIBINPUTS=$(srcdir): BSTINPUTS=$(srcdir):

if HAVE_LATEXMK

MBDating.pdf:
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/MBDating

else # !HAVE_LATEXMK

MBDating.pdf: \
    MBDating.tex \
    master_refs.bib \
    figures/ccby.pdf \
    figures/clockr.pdf \
    figures/cltree.pdf \
    figures/distribution.r \
    figures/fbdtree.pdf \
    figures/gamdir.pdf \
    figures/gamma.pdf \
    figures/nctree.pdf \
    figures/nodecali.pdf \
    figures/relaxcl.pdf \
    figures/treeage.pdf
	env $(envvars) $(PDFLATEX) MBDating
	env $(envvars) $(BIBTEX) MBDating
	env $(envvars) $(PDFLATEX) MBDating
	env $(envvars) $(PDFLATEX) MBDating

endif
else # !HAVE_LATEX

MBDating.pdf:
	echo 'Can not build the tutorial PDF since LaTeX was not detected' >&2
	exit 0

endif
