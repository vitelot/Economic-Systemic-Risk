include("extern.jl");
include("functions.jl");

function main(inputfile::String="data/test_list.csv", outputfile::String="data/output.csv")


    arrays_file = first(splitext(inputfile)) * "_arrays.jld2";
    
    # if file with arrays exists, load it with jld2
    if isfile(arrays_file)
        @info "Restoring working space from file \"$arrays_file\"";
        M,A = load(arrays_file, "M", "A");
        
    else 
        @info "Initializing the market according to input file \"$inputfile\"";
        M = initializeMarket(inputfile);
        @info "Building sparse adjacency matrices";
        A = buildArrays(M);
        @info "Saving matrices into file \"$arrays_file\"";
        jldsave(arrays_file; M,A);
    end

    @info "Calculating ESRI";
    esri = ESRI(M, A);

    @info "Saving the results into file \"$outputfile\"";
    saveESRI(M, esri, outputfile);

    return esri;
end

# esri = main("real_data/complete_edge_list.csv", "real_data/esri_complete.csv");
# esri = main("data/test_list.csv");
esri = main(ARGS...); # first arg is the input file, second is the output file
