# LaTeX and BibTex

## Install

In RedHat or CentOS System: <br />
```
yum install tetex
```

## Usage

Run latex on tex file to produce some output: <br />
```
latex doc.tex
bibtex doc.aux
latex doc.tex
latex doc.tex
```

View the output file (.dvi) with xdvi: <br />
```
xdvi doc.dvi &
```

To produce a PDF of this you simply run pdflatex instead of latex: <br />
```
pdflatex doc.tex
```

View the output file (.pdf) with evince: <br />
```
evince doc.pdf &
```

## Reference

[The Latex Project]

[The Latex Project]: https://www.latex-project.org
