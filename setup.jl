using CSV, DataFrames;

mutable struct Company
    nodeid::Int;
    # name::String;
    id::Int;
    sic::Int;
    nace::Int;
end

# graphfile    = "data/graph.jld2";
nodefile     = "data/nodes_ger.csv"; 
edgefile     = "data/edges_ger.csv"; 
sic2nacefile = "data/sic_to_nace2.csv";
essential_sector_file = "data/nace_essential.csv";

outfile      = "data/complete_edge_list.csv";

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