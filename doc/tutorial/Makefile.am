dist_doc_DATA=	hym.nex
dist_pdf_DATA=	MBDating.pdf

EXTRA_DIST=	src

CLEANFILES=	MBDating.aux \
		MBDating.bbl \
		MBDating.blg \
		MBDating.log \
		MBDating.out

DISTCLEANFILES=	MBDating.pdf

if BUILD_DOC
MBDating.pdf: \
    src/MBDating.tex \
    src/master_refs.bib \
    src/figures/ccby.pdf \
    src/figures/clockr.pdf \
    src/figures/cltree.pdf \
    src/figures/distribution.r \
    src/figures/fbdtree.pdf \
    src/figures/gamdir.pdf \
    src/figures/gamma.pdf \
    src/figures/nctree.pdf \
    src/figures/nodecali.pdf \
    src/figures/relaxcl.pdf \
    src/figures/treeage.pdf
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/MBDating
	BSTINPUTS=$(srcdir)/src BIBINPUTS=$(srcdir)/src $(BIBTEX) MBDating
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/MBDating
	TEXINPUTS=".:$(srcdir)/src:" $(PDFLATEX) $(srcdir)/src/MBDating
else
MBDating.pdf:
	echo 'Can not build the tutorial PDF since LaTeX was not detected' >&2
	exit 0
endif
