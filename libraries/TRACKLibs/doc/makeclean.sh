#!/bin/bash
name=$1
if [ -f "${name}.aux" ]; then rm ${name}.aux; fi
if [ -f "${name}.bbl" ]; then rm ${name}.bbl; fi
if [ -f "${name}.blg" ]; then rm ${name}.blg; fi
if [ -f "${name}.dvi" ]; then rm ${name}.dvi; fi
if [ -f "${name}.log" ]; then rm ${name}.log; fi
if [ -f "${name}.out" ]; then rm ${name}.out; fi
if [ -f "${name}.pdf" ]; then rm ${name}.pdf; fi
if [ -f "${name}.spl" ]; then rm ${name}.spl; fi
if [ -f "${name}.toc" ]; then rm ${name}.toc; fi
