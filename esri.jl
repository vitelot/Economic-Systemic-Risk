include("extern.jl");
include("functions.jl");

"""
    parseARGS(ARGS)

Parses command-line arguments using `ArgParse`.
Configuration options include input/output paths, max iterations (`--tmax`), shock scenarios (`--psi_mat`), and timeseries flags.
"""
function parseARGS(ARGS)
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--input", "-i"
            help = "The input file, formatted as a .csv with: supplier,customer,supplierNACE,customerNACE,weight,type"
            default = "data/test_list.csv" 
        "--output", "-o"
            help = "The desired output file"
            default = "data/output.csv"
        "--tmax", "-t"
            help = "Maximum iterations per scenario"
            arg_type = Int
            default = typemax(Int)
        "--psi_mat", "-p"
            help = "The psi_mat specifying the desired scenarios. A .csv with: scenario,firm,shocksize"
            default = ""
            arg_type = String
        "--timeseries"
            help = "Should the full timeseries of the ESRI calculation be returned?"
            action = :store_true
    end

    out = parse_args(ARGS, s)
    return out
end

"""
    main(ARGS)

Entry point for the script.
1. Parses arguments.
2. Initializes the `Market` (loading from CSV or cached `.jld2` arrays).
3. Constructs the shock scenarios (`psi_mat`).
4. Runs the `ESRI` calculation.
5. Saves the results to disk.
"""
function main(ARGS)
    ParsedARGS = parseARGS(ARGS)
    inputfile = ParsedARGS["input"]
    outputfile = ParsedARGS["output"]
    arrays_file = ""
    open(inputfile) do f
        arrays_file = bytes2hex(sha1(f)) * ".jld2"
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

    psi_mat = parsePsiMat(ParsedARGS["psi_mat"],M)

    @info "Calculating ESRI";
    esri = ESRI(M, A, psi_mat, ParsedARGS);

    @info "Saving the results into file \"$outputfile\"";
    saveESRI(M, esri, outputfile, ParsedARGS);

    return esri.esri;
end

esri = main(ARGS); # ARGS are parsed in the function. First 2 arguments should be input and output, then some options
