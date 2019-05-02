dist_doc_DATA=	hym.nex
dist_pdf_DATA=	MBDating.pdf

EXTRA_DIST=	src

DISTCLEANFILES=	MBDating.pdf

if HAVE_LATEX

CLEANFILES=	\
	MBDating.aux \
	MBDating.bbl \
	MBDating.blg \
	MBDating.fdb_latexmk \
	MBDating.fls \
	MBDating.log \
	MBDating.out

envvars=	\
	TEXINPUTS=$(srcdir)/src: \
	BIBINPUTS=$(srcdir)/src: \
	BSTINPUTS=$(srcdir)/src:

deps=	\
	$(srcdir)/src/MBDating.tex \
	$(srcdir)/src/figures/ccby.pdf \
	$(srcdir)/src/figures/clockr.pdf \
	$(srcdir)/src/figures/cltree.pdf \
	$(srcdir)/src/figures/fbdtree.pdf \
	$(srcdir)/src/figures/gamdir.pdf \
	$(srcdir)/src/figures/gamma.pdf \
	$(srcdir)/src/figures/nctree.pdf \
	$(srcdir)/src/figures/nodecali.pdf \
	$(srcdir)/src/figures/relaxcl.pdf \
	$(srcdir)/src/figures/treeage.pdf \
	$(srcdir)/src/master_refs.bib \
	$(srcdir)/src/sysbio.bst

if HAVE_LATEXMK

MBDating.pdf: $(deps)
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/src/MBDating

else # !HAVE_LATEXMK

MBDating.pdf: $(deps)
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/MBDating
	env $(envvars) $(BIBTEX) -terse MBDating
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/MBDating
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/MBDating

endif
else # !HAVE_LATEX

MBDating.pdf:
	echo 'Can not build the tutorial PDF since LaTeX was not detected' >&2
	exit 0

endif
