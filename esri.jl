include("extern.jl");
include("functions.jl");

function main(inputfile::String="data/input.csv", outputfile::String="data/output.csv")
    @info "Initializing the market according to input file \"$inputfile\"";
    M = initializeMarket(inputfile);

    @info "Building sparse adjacency matrices";
    A = buildArrays(M);

    nthreads = Threads.nthreads();
    if nthreads > 1
        @info "You set $nthreads cores to run the code in parallel. I'm impressed!";
        # @info "Calculating ESRI in parallel on $nthreads threads";
        esri = ESRI_parallel(M, A);
    else
        @info "Calculating ESRI sequentially";
        esri = ESRI(M, A);
    end

    @info "Saving the results into file \"$outputfile\"";
    saveESRI(M, esri, outputfile);

    return esri;
end

# esri = main("real_data/complete_edge_list.csv", "real_data/esri_complete.csv");
# esri = main("data/test_list.csv");
esri = main(ARGS...); # first arg is the input file, second is the output file
