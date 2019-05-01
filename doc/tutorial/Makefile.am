dist_doc_DATA=	hym.nex
dist_pdf_DATA=	MBDating.pdf

CLEANFILES=	\
	MBDating.aux \
	MBDating.bbl \
	MBDating.blg \
	MBDating.fdb_latexmk \
	MBDating.fls \
	MBDating.log \
	MBDating.out

DISTCLEANFILES=	MBDating.pdf

if HAVE_LATEX

envvars=	\
	TEXINPUTS=$(srcdir)/src: \
	BIBINPUTS=$(srcdir)/src: \
	BSTINPUTS=$(srcdir)/src:

if HAVE_LATEXMK

MBDating.pdf:
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/src/MBDating

else # !HAVE_LATEXMK

MBDating.pdf:
	env $(envvars) $(PDFLATEX) $(srcdir)/src/MBDating
	env $(envvars) $(BIBTEX) MBDating
	env $(envvars) $(PDFLATEX) $(srcdir)/src/MBDating
	env $(envvars) $(PDFLATEX) $(srcdir)/src/MBDating

endif
else # !HAVE_LATEX

MBDating.pdf:
	echo 'Can not build the tutorial PDF since LaTeX was not detected' >&2
	exit 0

endif
