function initializeMarket(file::String)::Market

    df = CSV.read(file, 
                types=Dict(:supplier => String, :customer => String), 
                DataFrame); 

    M = Market();

    C = M.Companies;
    Edges = M.Edges;
    NodeID = M.CompanyID;
    Sectors = M.Sectors;

    nodeid = 0;
    for r in eachrow(df)
        if !haskey(NodeID, r.supplier)
            nodeid += 1;
            NodeID[r.supplier] = nodeid;
        end
        if !haskey(NodeID, r.customer)
            nodeid += 1;
            NodeID[r.customer] = nodeid;
        end

        sid = NodeID[r.supplier];
        cid = NodeID[r.customer];
        snace = r.supplierNACE;
        cnace = r.customerNACE;
        w = r.weight;
        t = r.type;

        l = Edge(sid,cid,snace,cnace,w,t);

        get!(C, sid, Company(r.supplier, sid, snace, Edge[], Edge[], 0.0, 0.0));
        push!(C[sid].customers, l);

        get!(C, cid, Company(r.customer, cid, cnace, Edge[], Edge[], 0.0, 0.0));
        push!(C[cid].suppliers, l);

        push!(Edges, l);
    end

    for c in values(C)
        nace = c.nace;
        get!(Sectors, nace, Sector(nace));
        # push!(Sectors[nace].companies, c);
        # push!(Sectors[nace].company_ids, c.id);
        
        c.sout0 = sum([x.weight for x in c.customers]);
        c.sin0 = sum([x.weight for x in c.suppliers]);
    end
    for e in Edges
        t = e.type;
        cnace = e.customernace;
        snace = e.suppliernace;

        if t==2 # essential link
            push!(Sectors[cnace].essential, snace);
        elseif t==1 # non-essential link
            push!(Sectors[cnace].non_essential, snace);
        end

    end

    return M;
end

# function buildArrays(M::Market)::Arrays
#     # A = Arrays();

#     C = M.Companies;
#     Edges = M.Edges;

#     nrcomp = length(C);

#     # set up the weight matrix (fast way)
#     fromVect = [l.supplier for l in Edges];
#     toVect   = [l.customer for l in Edges];
#     weightVect = [l.weight for l in Edges];
#     W = sparse(fromVect, toVect, weightVect, nrcomp, nrcomp);

#     # W   = spzeros(nrcomp, nrcomp);
#     Λd1 = spzeros(nrcomp, nrcomp);
#     Λd2 = spzeros(nrcomp, nrcomp);
#     Λu  = spzeros(nrcomp, nrcomp);

#     sell = sparse(vec(sum(W, dims=2)));
#     # buy =  sparse(vec(sum(W, dims=1)));

#     for l in Edges
#         from = l.supplier;
#         to = l.customer;
#         val = W[from, to];
#         Λu[to, from] =  val / sell[from];
#         # println("$val $from $to $(buy[from]) $(sell[to])")
#     end

#     ### we build Λd1 and beta
#     β = spzeros(nrcomp);
#     # id = 2814;
#     for id in keys(C)
#         company = C[id];
#         D = Dict{Int,Vector{Float64}}(); # contains the supplier weights corresponding to the sector
#         F = Dict{Int,Vector{Int}}(); # contains the supplier id of the nodes corresponding to the sector
#         E = Dict{Int,Int}(); # link type i.e. if the supplier is essential
#         nonessential_weight_sum = 0.0;
#         for x in company.suppliers # cycle through the suppliers
#             suppliernace = x.suppliernace;
#             get!(D, suppliernace, Float64[]);
#             # append the weight to the corresponding vector addressed by the sector
#             push!(D[suppliernace], x.weight);
#             get!(F, suppliernace, Int[]);
#             # append the id to the corresponding vector addressed by the sector
#             push!(F[suppliernace], x.supplier);
#             get!(E, suppliernace, x.type);
#             nonessential_weight_sum += x.weight; # include ALL links
#             # x.type == 0 || (nonessential_weight_sum += x.weight); # exclude useless links
#             # x.type == 1 && (nonessential_weight_sum += x.weight); # include non-essential links only
#         end
#         all = beta = 0.0;
#         for sector in keys(D)
#             all += sum(D[sector]);
#             # if the sector is a necessary one for the node
#             if E[sector] == 2
#                 beta += sum(D[sector]);
#                 normalize!(D[sector],1);
#                 for i in eachindex(F[sector])
#                     Λd1[F[sector][i], id] = D[sector][i];
#                     # println("($(F[sector][i]), $nodeid) = $(D[sector][i])");
#                 end
#             elseif E[sector] == 1
#                 for i in eachindex(F[sector])
#                     Λd2[F[sector][i], id] = D[sector][i] / nonessential_weight_sum;
#                 end
#             end
#         end
#         if all>0
#             β[id] = beta/all;
#         end
#     end

#     # jldsave("data/arrays.jld2"; W, Λu, Λd1, Λd2, β)
    
#     return Arrays(Λu, Λd1, Λd2, β);
# end

"""
    buildArrays(M::Market)::Arrays

Constructs sparse matrices and vectors representing various relationships
within a market, based on company and edge data.

# Arguments
- `M::Market`: A Market object containing company and edge information.

# Returns
- `Arrays`: An object containing the constructed sparse matrices (Λu, Λd1, Λd2)
  and a sparse vector (β).
"""
function buildArrays(M::Market)::Arrays
    C = M.Companies
    Edges = M.Edges
    nrcomp = length(C)

    # Pre-allocate vectors for sparse matrix construction
    # W matrix
    W_I = Vector{Int}(undef, length(Edges))
    W_J = Vector{Int}(undef, length(Edges))
    W_V = Vector{Float64}(undef, length(Edges))

    # Λu matrix
    Λu_I = Vector{Int}(undef, length(Edges))
    Λu_J = Vector{Int}(undef, length(Edges))
    Λu_V = Vector{Float64}(undef, length(Edges))

    # Λd1 matrix
    Λd1_I = Int[]
    Λd1_J = Int[]
    Λd1_V = Float64[]

    # Λd2 matrix
    Λd2_I = Int[]
    Λd2_J = Int[]
    Λd2_V = Float64[]

    # β vector
    β_I = Int[]
    β_V = Float64[]

    # Populate W_I, W_J, W_V
    for (idx, l) in enumerate(Edges)
        W_I[idx] = l.supplier
        W_J[idx] = l.customer
        W_V[idx] = l.weight
    end
    W = sparse(W_I, W_J, W_V, nrcomp, nrcomp)

    # Calculate row sums for `sell` efficiently
    sell = vec(sum(W, dims=2))

    # Populate Λu_I, Λu_J, Λu_V using the calculated `sell` values
    # Ensure `sell[from]` is not zero to avoid division by zero.
    # If sell[from] is zero, it means the company has no outgoing edges,
    # so Λu[to, from] will effectively be zero.
    for (idx, l) in enumerate(Edges)
        from_node = l.supplier
        to_node = l.customer
        val = W_V[idx] # Use the original weight value from W_V
        if sell[from_node] != 0
            Λu_I[idx] = to_node
            Λu_J[idx] = from_node
            Λu_V[idx] = val / sell[from_node]
        else
            # If sell[from_node] is 0, this entry would be 0 or undefined,
            # so we can skip adding it to the sparse list or set its value to 0.
            # Here, we set to 0, but sparse matrix construction will omit 0s.
            Λu_I[idx] = to_node
            Λu_J[idx] = from_node
            Λu_V[idx] = 0.0
        end
    end
    # Filter out zero values before constructing the sparse matrix for efficiency
    # Although sparse() itself handles zeros, explicitly filtering might be marginally faster
    # if many zeros are produced. Given the nature of the division, many won't be zero.
    # Keep as is, let sparse() handle zeros during construction.
    Λu = sparse(Λu_I, Λu_J, Λu_V, nrcomp, nrcomp)

    # Process companies for Λd1, Λd2, and β
    # We iterate through company IDs to ensure all companies are considered,
    # even if they have no suppliers (though current logic would yield no entries for them).
    for id in keys(C)
        company = C[id]
        
        # Use mutable dictionaries for collecting supplier data by sector
        D_sector_weights = Dict{Int, Vector{Float64}}() # supplier weights by sector
        F_sector_ids = Dict{Int, Vector{Int}}()         # supplier IDs by sector
        E_sector_type = Dict{Int, Int}()                 # sector type (essentiality)

        nonessential_weight_sum = 0.0

        for x in company.suppliers
            supplier_nace = x.suppliernace
            
            get!(D_sector_weights, supplier_nace, Float64[])
            push!(D_sector_weights[supplier_nace], x.weight)

            get!(F_sector_ids, supplier_nace, Int[])
            push!(F_sector_ids[supplier_nace], x.supplier)
            
            E_sector_type[supplier_nace] = x.type # Overwrites if multiple suppliers in same sector, assuming type is consistent

            # Summing ALL weights for `nonessential_weight_sum` based on the original logic
            nonessential_weight_sum += x.weight
        end

        total_weight_all_suppliers = 0.0
        total_weight_essential_suppliers = 0.0

        for sector in keys(D_sector_weights)
            current_sector_weights = D_sector_weights[sector]
            current_sector_supplier_ids = F_sector_ids[sector]
            sector_type = get(E_sector_type, sector, 0) # Default to 0 if not found, though should always be there

            sector_sum_weight = sum(current_sector_weights)
            total_weight_all_suppliers += sector_sum_weight

            if sector_type == 2 # Essential sector
                total_weight_essential_suppliers += sector_sum_weight
                
                # Normalize weights within this essential sector
                if(sector_sum_weight > 0.0)
                    normalized_weights = current_sector_weights ./ sector_sum_weight
                else
                    # there are all zeros in the vector since sector_sum_weight is zero
                    normalized_weights = current_sector_weights;
                end
                
                for i in eachindex(current_sector_supplier_ids)
                    push!(Λd1_I, current_sector_supplier_ids[i])
                    push!(Λd1_J, id)
                    push!(Λd1_V, normalized_weights[i])
                end
            elseif sector_type == 1 # Non-essential sector
                # Note: Original code used `nonessential_weight_sum` for normalization here.
                # If nonessential_weight_sum is 0, these values will be 0.
                if nonessential_weight_sum != 0
                    for i in eachindex(current_sector_supplier_ids)
                        push!(Λd2_I, current_sector_supplier_ids[i])
                        push!(Λd2_J, id)
                        push!(Λd2_V, current_sector_weights[i] / nonessential_weight_sum)
                    end
                end
            end
        end

        # Calculate β for the current company
        if total_weight_all_suppliers > 0
            push!(β_I, id)
            push!(β_V, total_weight_essential_suppliers / total_weight_all_suppliers)
        end
    end

    Λd1 = sparse(Λd1_I, Λd1_J, Λd1_V, nrcomp, nrcomp)
    Λd2 = sparse(Λd2_I, Λd2_J, Λd2_V, nrcomp, nrcomp)
    β = sparse(β_I, ones(Int, length(β_I)), β_V, nrcomp, 1) # β is a column vector

    return Arrays(Λu, Λd1, Λd2, β)
end

function marketShare(M::Market, Q::DynamicalQuantities)::Nothing
    C = M.Companies;
    # Sectors = M.Sectors;
    marketshare = Q.marketshare; #spzeros(nrcomp);
    hd = Q.hd;
    
    volumesector = Dict{Int,Float64}();
    for c in values(C)
        volumesector[c.nace] = get(volumesector, c.nace, 0.0) + c.sout0 * hd[c.id];
    end
    for company in values(C)
        sout0 = company.sout0;
        vol_sec = volumesector[company.nace];
        
        if sout0 > 0
            if vol_sec > 0.0
                marketshare[company.id] = min(1.0, sout0 / vol_sec);
            else
                marketshare[company.id] = 1.0;
            end
        else
            marketshare[company.id] = 0.0;
        end
    end
    return;
end

function upStream(company::Company, A::Arrays, hu::Vector{Float64})::Float64
    company.sout0 == 0.0 && return 1.0; # no customers -> no upstream shock
    D_u = 0.0;
    id = company.id;
    for e in company.customers
        cid = e.customer;
        D_u += A.lambda_u[cid, id] * hu[cid];
    end
    return D_u;
end

function downStream(company::Company, A::Arrays, Q::DynamicalQuantities)::Tuple{Float64,Float64}

    marketshare = Q.marketshare;
    hd = Q.hd;

    id = company.id;

    D = Dict{Int,Float64}(); # partial sums by essential sector, i.e., Π_ik in the paper
    D_ne = 0.0; # non-essential contribution 

    for e in company.suppliers
        t = e.type;
        snace = e.suppliernace;
        sid = e.supplier;
        if t == 2
            D[snace] = get(D, snace, 0.0) + marketshare[sid] * A.lambda_d1[sid, id] * (1.0 - hd[sid]);
        elseif t==1
            D_ne += marketshare[sid] * A.lambda_d2[sid, id] * (1.0 - hd[sid]);
        end
    end
    # println(id)
    essentials = 1.0 - maximum(values(D), init=0.0);
    non_essentials = 1.0 - D_ne;

    return essentials, non_essentials;
end

function oneStep(M::Market, A::Arrays, Q::DynamicalQuantities)::Float64
    C = M.Companies;
    # Sectors = M.Sectors;
    hd = Q.hd;
    hu = Q.hu;
    ψ = Q.psi;

    newhd = Q.newhd ; #copy(hd); # use similar later
    newhu = Q.newhu; #copy(hu);

    marketShare(M,Q);

    # company = C[1302];
    for company in values(C)
        id = company.id;
        essentials, non_essentials = downStream(company, A, Q);
        newhd[id] = minimum((essentials, non_essentials, ψ[id]));
        
        D_u = upStream(company, A, hu);
        newhu[id] = min(D_u, ψ[id]);
        # println("$id, $essentials, $non_essentials, $D_u, $(newhd[id]), $(newhu[id])");
    end

    error = max( maximum( abs.(hd .- newhd) ), maximum( abs.(hu .- newhu) ) );

    # garbage collector friendly: copy vectors without changing Q
    Q.hd .= newhd;
    Q.hu .= newhu;

    return error;
end

function ESRI(M::Market, A::Arrays)::Vector{Float64}
    
    nrcomp = length(M.Companies);
    
    @info "Initializing dynamical quantities for $nrcomp companies";
    Q = DynamicalQuantities(nrcomp);

    esri = Vector{Float64}(undef, nrcomp);
    # h    = Vector{Float64}(undef, nrcomp);
    total_volume = sum([x.sout0 for x in values(M.Companies)]);
    u = ones(nrcomp);
    # @showprogress dt=1 desc="Computing..." for t in 1:tmax
        @showprogress for i in sort(collect(keys(M.Companies)))
            Q.psi .= u; Q.psi[i] = 0.0;
            Q.hd .= u; Q.hu .= u;
            Q.newhd .= u; Q.newhu .= u;

            err = 1.0;
            while err > 1e-2
                # println("------------------\n$i $err");
                err = oneStep(M,A,Q);
            end
            h = 1.0 .- min.(Q.hd, Q.hu);
            esri[i] = sum([x.sout0 * h[x.id] for x in values(M.Companies)]) / total_volume;
            # println("ESRI[$(M.Companies[i].name)] = $(esri[i])");
            @assert !isnan(esri[i]);
        end
    return esri;
end

function ESRI_parallel(M::Market, A::Arrays)::Vector{Float64}
    
    nthreads = Threads.nthreads();
    if nthreads == 1
        println("You have only one thread selected. If you have more, pls run julia with --threads n");
        exit();
    end
    
    Results = DataFrame(index=Int[], esri=Float64[]);
    VR = Vector{DataFrame}(undef, nthreads);
    VQ = Vector{DynamicalQuantities}(undef, nthreads);
    for i = 1:nthreads
        VR[i] = copy(Results);
        VQ[i] = DynamicalQuantities(length(M.Companies));
    end
    # @info 1
    nrcomp = length(M.Companies);
    total_volume = sum([x.sout0 for x in values(M.Companies)]);
    u = ones(nrcomp);

    Threads.@threads for i in collect(keys(M.Companies))
        tid = Threads.threadid();
        VQ[tid].psi .= u; VQ[tid].psi[i] = 0.0;
        VQ[tid].hd .= u; VQ[tid].hu .= u;
        err = 1.0;
        while err > 1e-2
            err = oneStep(M,A,VQ[tid]);
        end
        # println("Calculating esri for firm \"$(M.Companies[i].name)\" on thread $tid ");
        h = 1.0 .- min.(VQ[tid].hd, VQ[tid].hu);
        esri = sum([x.sout0 * h[x.id] for x in values(M.Companies)]) / total_volume;
        push!(VR[tid], (i, esri));
    end

    df = sort(vcat(VR...), :index);
    return df.esri;
end

function saveESRI(M::Market, esri::Vector{Float64}, outfile::String)
    companies = [x.name for x in sort(collect(values(M.Companies)), by=x->x.id)];
    dfout = DataFrame(company=companies, esri=esri);
    CSV.write(outfile, dfout);
end