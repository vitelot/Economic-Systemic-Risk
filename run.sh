#! /bin/sh

echo "Enter the number of threads you wish to use:"
read nthreads
echo 

if [ -z "$nthreads" ]; then
    nthreads=1
fi

# Define defaults for clarity
default_in="data/test_list.csv"
default_out="data/output.csv"

if [ $# -eq 0 ]; then
    echo "No filenames provided. Continuing with default $default_in and $default_out"
    in="$default_in"
    out="$default_out"
    psi=""
elif [ $# -eq 1 ]; then
    echo "No output filename provided. Continuing with default $default_out"
    in="$1"
    out="$default_out"
    psi=""
elif [ $# -eq 2 ]; then
    echo "Running ESRI in Julia with input file ${1} and output ${2}"
    in="$1"
    out="$2"
    psi=""
elif [ $# -eq 3 ]; then
    echo "Running ESRI in Julia with input file ${1}, output ${2}, and psi mat ${3}"
    in="$1"
    out="$2"
    psi="$3"
fi


echo "julia --project --threads $nthreads main.jl -i $in -o $out -p \"$psi\""
julia --project --threads $nthreads main.jl -i "$in" -o "$out" -p "$psi"