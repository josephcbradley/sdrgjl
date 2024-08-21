#!/bin/bash -l
julia --threads=16 --color=yes --startup-file=no --history-file=no --project=. src/full_run.jl
cd tex 
pdflatex --file-line-error --shell-escape main.tex
bibtex main.aux
pdflatex --file-line-error --shell-escape main.tex
cd ..
