#!/bin/bash

rm main.aux  main.toc main.snm   main.nav  main.vrb   main.out   main.aux   main.pdf  main.log main.blg main.bbl

pdflatex -shell-escape main
bibtex main
pdflatex -shell-escape main
pdflatex -shell-escape main
#dvips main.dvi 
#ps2pdf main.ps
#okular main.pdf

#pdflatex main.tex
