using CSV, DataFrames;

"""
    Company

A temporary mutable struct used during the data preprocessing phase to map firm attributes before generating the final edge list.

# Fields
- `nodeid::Int`: Sequential internal index assigned during reading.
- `id::Int`: Original identifier from the raw data.
- `sic::Int`: Standard Industrial Classification (SIC) code (raw input).
- `nace::Int`: NACE sector code (mapped from SIC).
"""
mutable struct Company
    nodeid::Int;
    # name::String;
    id::Int;
    sic::Int;
    nace::Int;
end

# graphfile    = "data/graph.jld2";
nodefile     = "real_data/nodes_ger.csv"; 
edgefile     = "real_data/edges_ger.csv"; 
sic2nacefile = "real_data/sic_to_nace2.csv";
essential_sector_file = "real_data/nace_essential.csv";

outfile      = "real_data/complete_edge_list.csv";

"""
    main()

The primary preprocessing routine. It transforms raw real-world data files into the standardized edge list format required by `esri.jl`.

# Workflow
1. **Load Data**: Reads raw node attributes, edge lists, SIC-to-NACE conversion tables, and sector essentiality matrices.
2. **Map Sectors**: Converts US SIC codes to European NACE codes for every firm.
3. **Classify Links**: Iterates through the raw edge list and assigns a `type` to each connection based on the sector-to-sector relationship (defined in `nace_essential.csv`).
4. **Export**: Saves the fully processed network (Supplier, Customer, NACEs, Weight, Type) to `real_data/complete_edge_list.csv`.
"""
function main()
    @info "Loading node file";
    dfnode = CSV.read(nodefile, select=[:CompanyID, :SIC_PrimaryIndustryCode], DataFrame);
    dfnode[!,2] = coalesce.(dfnode[!,2], 9999);
    # dfnode  = CSV.read(nodefile, DataFrame)

    @info "Loading edge file";
    dfedge  = CSV.read(edgefile, DataFrame)

    @info "Loading SIC to NACE2 conversion table";
    df_s2n = CSV.read(sic2nacefile, DataFrame);
    D_s2n = Dict(Pair.(df_s2n.SIC, df_s2n.NACE2));

    @info "Loading sector relationships";
    dfrel = CSV.read(essential_sector_file, DataFrame);
    S = Dict{Tuple{Int,Int},Int}();
    for r in eachrow(dfrel)
        S[(r[1],r[2])] = r[3];
    end

    C = Dict{Int, Company}();

    nodeid = 0;
    for r in eachrow(dfnode)
        nodeid += 1;
        id = r.CompanyID;
        sic = r.SIC_PrimaryIndustryCode;
        nace = D_s2n[sic÷100]; # take the first 2 digits)
        c = Company(nodeid, id, sic, nace);
        C[id] = c;
    end 

    dfout = DataFrame(supplier=Int[], customer=Int[], 
                    # supplierID=Int[], customerID=Int[],
                    supplierNACE=Int[], customerNACE=Int[],
                    weight=Float64[], type=Int[]);

    # since there are companies that do not show up in the edge list we 
    # assign their ID from the edge list, not from the node list
    
    for r in eachrow(dfedge)
        suppliernace = C[r.SuppID].nace;
        customernace = C[r.CustID].nace;
        # suppliernace==customernace || 
            push!(dfout, (r.SuppID, r.CustID, suppliernace, customernace, 1.0, S[(suppliernace,customernace)]));
    end
    @info "Saving the complete edge list to file \"$outfile\"";
    CSV.write(outfile, dfout);
end

main();