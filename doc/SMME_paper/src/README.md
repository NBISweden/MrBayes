# Files for building "Bayesian Analysis of Molecular Evolution using MrBayes"

- Last modified: fre jul 26, 2019  07:04
- Sign: Johan Nylander


## Description

* The files in this directory are "reversed engineered" to be able to build
the LaTeX file `Huelsenbeck_Ronquist_SMME_2005.tex`. The build process works, 
but see **Issues** below.

* This is a **draft version** to what was published as

> Huelsenbeck J.P., Ronquist F. (2005) Bayesian Analysis of Molecular Evolution
> Using MrBayes. Pp 183-226 in: Statistical Methods in Molecular Evolution. Statistics
> for Biology and Health. Springer, New York, NY.
>
> DOI: <https://doi.org/10.1007/0-387-27733-1_7>
>
> URL: <https://link.springer.com/chapter/10.1007/0-387-27733-1_7>

* The style class file `svmult.cls` is from Springer, and was downloaded on 19 July 2019 from

[https://sv.sharelatex.com/templates/books/springer's-edited-book-svmult](https://sv.sharelatex.com/templates/books/springer's-edited-book-svmult)


## Files included 

    Huelsenbeck_Ronquist_SMME_2005.tex
    Huelsenbeck_Ronquist_SMME_2005.bib
    svmult.cls
    fig1.png  fig2.png  fig3.png  fig4.png
    fig5.png  fig6.png  fig7.png  fig8.png
    fig9.png  fig10.png fig11.png fig12.png
    README.md

## Compile

    pdflatex Huelsenbeck_Ronquist_SMME_2005.tex
    bibtex Huelsenbeck_Ronquist_SMME_2005
    pdflatex Huelsenbeck_Ronquist_SMME_2005

Or, if `latexmk` is available:

    latexmk -pdf Huelsenbeck_Ronquist_SMME_2005.tex


## Issues

- The bibliography list is not yet complete. See the file
  `Huelsenbeck_Ronquist_SMME_2005.bib`.

- Figures are not original versions.
  I did not have access to the original image material, so all PNG images
  are screenshots(!) from the original publication.

- "`Overfull \hbox`" issues -- several.
  I changed font size for some of the verbatim listings in the Appendix 2 to
  accommodate some warnings.

- "`Underfull \hbox`" issues -- several.
  No action taken.

- `LaTeX Font Warning: Font shape `OT1/lmr/bx/sc' undefined`.
  No action taken.

- Use of alternative command:
```
        Package amsmath Warning: Foreign command \over;
        (amsmath)                \frac or \genfrac should be used instead
        (amsmath)                 on input line 108.
```
  This can be adressed by using, for example, `\frac{a}{b}` instead of `{a \over b}`.
  No action taken.
    
