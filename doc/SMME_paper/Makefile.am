dist_pdf_DATA=	Huelsenbeck_Ronquist_SMME_2005.pdf

EXTRA_DIST=	src

DISTCLEANFILES=	Huelsenbeck_Ronquist_SMME_2005.pdf

CLEANFILES=	\
	Huelsenbeck_Ronquist_SMME_2005.aux \
	Huelsenbeck_Ronquist_SMME_2005.bbl \
	Huelsenbeck_Ronquist_SMME_2005.blg \
	Huelsenbeck_Ronquist_SMME_2005.fdb_latexmk \
	Huelsenbeck_Ronquist_SMME_2005.fls \
	Huelsenbeck_Ronquist_SMME_2005.log \
	Huelsenbeck_Ronquist_SMME_2005.idx \
	Huelsenbeck_Ronquist_SMME_2005.ilg \
	Huelsenbeck_Ronquist_SMME_2005.ind \
	Huelsenbeck_Ronquist_SMME_2005.toc

envvars=	\
	TEXINPUTS=$(srcdir)/src: \
	BIBINPUTS=$(srcdir)/src: \
	BSTINPUTS=$(srcdir)/src:

deps=	\
	$(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005.bib \
	$(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005.tex \
	$(srcdir)/src/fig1.png \
	$(srcdir)/src/fig2.png \
	$(srcdir)/src/fig3.png \
	$(srcdir)/src/fig4.png \
	$(srcdir)/src/fig5.png \
	$(srcdir)/src/fig6.png \
	$(srcdir)/src/fig7.png \
	$(srcdir)/src/fig8.png \
	$(srcdir)/src/fig9.png \
	$(srcdir)/src/fig10.png \
	$(srcdir)/src/fig11.png \
	$(srcdir)/src/fig12.png \
	$(srcdir)/src/svmult.cls

if HAVE_LATEX
if HAVE_LATEXMK

pdf: Huelsenbeck_Ronquist_SMME_2005.pdf

Huelsenbeck_Ronquist_SMME_2005.pdf: $(deps)
	env $(envvars) $(LATEXMK) -quiet -pdf $(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005

else # !HAVE_LATEXMK

# The order of pdflatex and bibtex runs below (and the number of runs)
# mirrors what latexmk seems to want to do.
Huelsenbeck_Ronquist_SMME_2005.pdf: $(deps)
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005
	env $(envvars) $(BIBTEX) -terse Huelsenbeck_Ronquist_SMME_2005
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005
	env $(envvars) $(BIBTEX) -terse Huelsenbeck_Ronquist_SMME_2005
	env $(envvars) $(PDFLATEX) --interaction=batchmode $(srcdir)/src/Huelsenbeck_Ronquist_SMME_2005

endif
else # !HAVE_LATEX

Huelsenbeck_Ronquist_SMME_2005.pdf:
	@echo 'Can not build the Huelsenbeck and Ronquist 2005 PDF since LaTeX was not detected' >&2
	@exit 0

endif
