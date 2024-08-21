#!/bin/bash -l

julia --threads=8 --color=yes --startup-file=no --history-file=no --project=. src/part_run.jl