# -*- mode: makefile -*-
#========================================
# Makefile for a static HTML Site
#
# Author: Frédéric Haziza <daz@bils.se>
#========================================

CONTENT_FILES := $(wildcard *.html)
GENERATED_FILES := $(CONTENT_FILES:%.html=../%.html)

all: $(GENERATED_FILES)

../%.html: %.html Makefile
	@echo "Creating the $* page [in $@]"
	@echo '<!DOCTYPE html>' > $@
	@echo '<html lang=en><head><meta charset="utf-8"><title>MrBayes | $*</title>' >> $@
	@echo '<link rel=stylesheet href="css/style.css"/>' >> $@
	@echo '<link rel=stylesheet href="css/handheld.css" media="only screen and (min-device-width:320px) and (max-device-width:568px)" />' >> $@
	@echo '<meta name="viewport" content="width=device-width, initial-scale=1.0" />' >> $@
	@echo '</head><body>' >> $@
	@echo '<header><h1>MrBayes: Bayesian Inference of Phylogeny</h1>' >> $@
	@echo '<nav>' >> $@
	@echo '<a href="./index.html"' >> $@
	@if [ "$*" = "index" ]; then echo 'class=selected' >> $@;fi
	@echo '>Home</a>' >> $@
	@echo '<a href="./download.html"' >> $@
	@if [ "$*" = "download" ]; then echo 'class=selected' >> $@;fi
	@echo '>Download</a>' >> $@
	@echo '<a href="./manual.html"' >> $@
	@if [ "$*" = "manual" ]; then echo 'class=selected' >> $@;fi
	@echo '>Manual</a>' >> $@
	@echo '<a href="./issues.html"' >> $@
	@if [ "$*" = "issues" ]; then echo 'class=selected' >> $@;fi
	@echo '>Bug Report</a>' >> $@
	@echo '<a href="./authors.html"' >> $@
	@if [ "$*" = "authors" ]; then echo 'class=selected' >> $@;fi
	@echo '>Authors</a>' >> $@
	@echo '<a href="./links.html"' >> $@
	@if [ "$*" = "links" ]; then echo 'class=selected' >> $@;fi
	@echo '>Links</a>' >> $@
	@echo '</nav>' >> $@
	@echo '<aside>' >> $@
	@echo '<a href="http://www.vr.se/"><img src="./img/vr.png" alt="Vetenskapsradet"></a>' >> $@
	@echo '<a href="http://www.nih.gov/"><img src="./img/nih.gif" alt="National Institutes of Health"></a>' >> $@
	@echo '<a href="http://www.nsf.gov/"><img src="./img/nsf.gif" alt="National Science Foundation"></a>' >> $@
	@echo '<a href="http://www.nrm.se/"><img src="./img/nrm.gif" alt="Naturhistoriska Riksmuseet"></a>' >> $@
	@echo '<a href="http://nbis.se/"><img src="./img/nbis.svg" alt="National Bioinformatics Infrastructure Sweden (NBIS)"></a>' >> $@
	@echo '</aside>' >> $@
	@echo '</header>' >> $@
	@echo '<main id=$*>' >> $@
	@cat $< >> $@
	@echo '</main>' >> $@
	@echo '</body></html>' >> $@
	@tidy -config tidy.conf $@

echo:
	@echo $(GENERATED_FILES)
clean:
	@echo "Cleaning the generated HTML files"
	@rm -rf $(GENERATED_FILES)
