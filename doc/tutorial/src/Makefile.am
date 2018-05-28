all:	MBDating.pdf

EXTRA_DIST=	MBDating.tex \
		master_refs.bib \
		figures \
		sysbio.bst

CLEANFILES=	MBDating.aux \
		MBDating.bbl \
		MBDating.blg \
		MBDating.log \
		MBDating.out

DISTCLEANFILES=	MBDating.pdf

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
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/MBDating
	BSTINPUTS=$(srcdir) BIBINPUTS=$(srcdir) $(BIBTEX) MBDating
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/MBDating
	TEXINPUTS=".:$(srcdir):" $(PDFLATEX) $(srcdir)/MBDating
