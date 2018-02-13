#!/bin/bash
name=$1
latex $name
bibtex $name
latex $name
latex $name
pdflatex $name
