#! /bin/sh

if [ $# = 0 ]; then
    echo "No filenames provided. Continuing with default data/input.csv and data/output.csv"
elif [ $# = 1 ]; then
    echo "No output filename provided. Continuing with default data/output.csv"
elif [ $# = 2 ]; then
    echo "Running ESRI in Julia with input file ${1} and output ${2}"
fi

echo "julia --project esri.jl $@"
julia --project esri.jl $@