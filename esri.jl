include("extern.jl");
include("functions.jl");

using ArgParse, SHA

function parseARGS(ARGS)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input"
            help = "The input file, formatted as a .csv with: supplier,customer,supplierNACE,customerNACE,weight,type"
            default = "data/test_list.csv" 
        "--output"
            help = "The desired output file"
            default = "data/output.csv"
        "--tmax", "-t"
            help = "Maximum iterations per scenario"
            arg_type = Int
            default = typemax(Int)
        "--psi_mat", "-p"
            help = "The psi_mat specifying the desired scenarios. A .csv with: scenario,firm,shocksize"
            default = 0
        "--timeseries"
            help = "Should the full timeseries of the ESRI calculation be returned? Not yet implemented"
            action = :store_true
    end

    out = parse_args(ARGS, s)
    if out["timeseries"]
        error("Not implemented: timeseries")
    end
    if out["psi_mat"]!=0
        error("Not implemented: psi_mat")
    end
    return out
end

function main(ARGS)
    ParsedARGS = parseARGS(ARGS)
    inputfile = ParsedARGS["input"]
    outputfile = ParsedARGS["output"]
    arrays_file = ""
    open(inputfile) do f
        arrays_file = bytes2hex(sha512(f)) * ".jld2"
    end
    
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
    esri = ESRI(M, A, ParsedARGS);

    @info "Saving the results into file \"$outputfile\"";
    saveESRI(M, esri, outputfile);

    return esri;
end

# esri = main("real_data/complete_edge_list.csv", "real_data/esri_complete.csv");
# esri = main("data/test_list.csv");
esri = main(ARGS); # first arg is the input file, second is the output file
