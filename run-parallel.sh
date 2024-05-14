#! /bin/sh

echo "Enter the number of threads you wish to use:"
read nthreads
echo 

if [ $# = 0 ]; then
    echo "No filenames provided. Continuing with default data/input.csv and data/output.csv"
elif [ $# = 1 ]; then
    echo "No output filename provided. Continuing with default data/output.csv"
elif [ $# = 2 ]; then
    echo "Running ESRI in Julia with input file ${1} and output ${2}"
fi

echo "julia --project --threads $nthreads esri.jl $@"
julia --project --threads $nthreads esri.jl $@