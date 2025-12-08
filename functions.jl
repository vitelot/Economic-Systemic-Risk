"""
    initializeMarket(file::String)::Market

Reads the input CSV file (edge list) and constructs the `Market` object.
It initializes `Company` nodes, `Edge` connections, and `Sector` classifications,
identifying essential and non-essential links based on the input columns.
"""
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

function initSectorVolumes!(M::Market, Q::DynamicalQuantities)
    for c in values(M.Companies)
        Q.initial_sector_volumes[c.nace] = get(Q.initial_sector_volumes, c.nace, 0.0) + c.sout0;
    end

end

"""
    parsePsiMat(file::String, M::Market)

Parses the scenario matrix file (`psi_mat`) defining the shocks.
- If `file` is empty, it returns a diagonal sparse identity matrix (implying single-firm shock scenarios for every firm).
- Otherwise, it reads the CSV to construct a sparse matrix where columns represent scenarios and rows represent firms.
"""
function parsePsiMat(file::String, M::Market)
    nfirms = length(M.Companies)
    if file==""
        @info "Using diagonal psi_mat"
        return sparse(I(nfirms))
    end
    @info "Reading provided psi_mat at $file"
    data = CSV.read(file,DataFrame,types=Dict(1=>Int,2=>String,3=>Float64))
    firmids = getindex.(Ref(M.CompanyID),strip.(data[:,2])) # strip removes leading and trailing whitespaces. In case the csv is a bit broken
    nscenarios = maximum(data[:,1])
    psi_mat = sparse(firmids,data[:,1],data[:,3],nfirms,nscenarios) 
    return psi_mat
end

"""
    buildArrays(M::Market)::Arrays

Converts the object-oriented `Market` graph into the sparse algebraic structures required for the linear algebra operations.

# Returns
`Arrays` struct containing:
- `lambda_u`: Upstream adjacency matrix (customer to supplier).
- `lambda_d1`: Downstream adjacency matrix for **essential** inputs.
- `lambda_d2`: Downstream adjacency matrix for **non-essential** inputs.
- `beta`: Vector representing the share of essential inputs for each firm.
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

"""
    marketShare(M::Market, Q::DynamicalQuantities)::Nothing

Updates the dynamic market share of each company within its sector.
Calculated as the ratio of the company's current output (s_{out,0} ⋅ h_d) 
to the total current output of its sector.
"""
function marketShare(M::Market, Q::DynamicalQuantities)::Nothing
    C = M.Companies;
    # Sectors = M.Sectors;
    marketshare = Q.marketshare; #spzeros(nrcomp);
    hd = Q.hd;
    
    sector_volumes = copy(Q.initial_sector_volumes);
    # for c in values(C)
    #     sector_volumes[c.nace] = get(sector_volumes, c.nace, 0.0) + c.sout0 * hd[c.id];
    # end
    for id in Q.changed_firms
        c = M.Companies[id];
        sector_volumes[c.nace] += c.sout0 * (hd[c.id] - 1.0);
    end

    for company in values(C)
        sout0 = company.sout0;
        vol_sec = sector_volumes[company.nace];
        cid = company.id;
        if sout0 > 0
            if vol_sec > 0.0
                marketshare[cid] = min(1.0, sout0 / vol_sec);
            else
                marketshare[cid] = 1.0;
            end
        else
            marketshare[cid] = 0.0;
        end
        # println("$cid $(marketshare[cid])");
    end
    return;
end

"""
    upStream(company::Company, A::Arrays, hu::Vector{Float64})::Float64

Calculates the upstream demand shock (D_u) for a specific `company`.
This aggregates the health of customers (h_u) weighted by the upstream matrix Λ_u.
"""
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

"""
    downStream(company::Company, A::Arrays, Q::DynamicalQuantities)::Tuple{Float64,Float64}

Calculates the downstream supply availability for a specific `company`.

# Returns
A tuple `(essentials, non_essentials)`:
- `essentials`: The available supply from essential inputs (using Λ_d1).
- `non_essentials`: The available supply from non-essential inputs (using Λ_d2).
"""
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

"""
    oneStep(M::Market, A::Arrays, Q::DynamicalQuantities)::Float64

Performs a single iteration of the fixed-point algorithm to update production levels.
1. Updates market shares.
2. Computes new downstream (h_d) and upstream (h_u) levels for all firms.
3. Updates `Q` in place.

# Returns
The maximum error (Chebyshev distance) between the previous and current state, used for convergence checking.
"""
function oneStep(M::Market, A::Arrays, Q::DynamicalQuantities)::Float64
    C = M.Companies;
    # Sectors = M.Sectors;
    hd = Q.hd;
    hu = Q.hu;
    ψ = Q.psi;

    newhd = Q.newhd ;
    newhu = Q.newhu;

    marketShare(M,Q);

    error = 0.0;
    for company in values(C)
        id = company.id;
        essentials, non_essentials = downStream(company, A, Q);
        newhd[id] = min(essentials, non_essentials, ψ[id]);
        
        D_u = upStream(company, A, hu);
        newhu[id] = min(D_u, ψ[id]);
        # println("$id, $essentials, $non_essentials, $D_u, $(newhd[id]), $(newhu[id])");
        c_error = max(error, abs(newhd[id]-hd[id]), abs(newhu[id]-hu[id]));
        if c_error > 1e-9 # firm's status changed
            push!(Q.changed_firms, id);
            # println("Changed: $id");
        end
        error = c_error;
    end

    # error = max( maximum( abs.(hd .- newhd) ), maximum( abs.(hu .- newhu) ) );

    # garbage collector friendly: copy vectors without changing Q
    Q.hd .= newhd;
    Q.hu .= newhu;

    return error;
end

"""
    ESRI(M::Market, A::Arrays, psi_mat::SparseMatrixCSC, ParsedARGS)::DataFrame

The main simulation engine.
It iterates through scenarios defined in `psi_mat` (columns) using multi-threading (`Threads.@threads`).
For each scenario, it converges the `oneStep` function until the error is below a threshold or `tmax` is reached.

# Returns
A `DataFrame` containing the index of the scenario, the calculated ESRI value, and the number of iterations required.
"""
function ESRI(M::Market, A::Arrays, psi_mat::SparseMatrixCSC, ParsedARGS)::DataFrame
    tmax = ParsedARGS["tmax"] 
 
    # FIX: Use maxthreadid() if available to handle non-contiguous thread IDs
    # This prevents the "BoundsError at index [3]" when nthreads is 2
    max_tid = isdefined(Threads, :maxthreadid) ? Threads.maxthreadid() : Threads.nthreads()
    
    @info "Using $max_tid threads"

    Results = DataFrame(index=Int[], esri=Float64[], t = Int[]);
    
    # Allocate vectors up to the highest possible Thread ID
    VR = Vector{DataFrame}(undef, max_tid);
    VQ = Vector{DynamicalQuantities}(undef, max_tid);
    
    # Initialize all potential slots
    for i = 1:max_tid
        VR[i] = copy(Results);
        VQ[i] = DynamicalQuantities(length(M.Companies));
        initSectorVolumes!(M, VQ[i]);
    end
    nrcomp = length(M.Companies);
    total_volume = sum([x.sout0 for x in values(M.Companies)]);
    u = ones(nrcomp);

    Threads.@threads for i in axes(psi_mat,2) # go through all firms unless a psi scenario loaded
        t = 0
        tid = Threads.threadid();
        indices = nzrange(psi_mat,i)
        firms = rowvals(psi_mat)[indices]
        psi = 1 .- nonzeros(psi_mat)[indices]
        
        # println("Calculating esri for firm \"$(M.Companies[i].name)\" on thread $tid ");
        VQ[tid].psi .= u; VQ[tid].psi[firms] = psi;
        VQ[tid].hd .= u; VQ[tid].hu .= u;

        empty!(VQ[tid].changed_firms);
        # pprintln(VQ[tid]);

        err = 1.0;
        while (err > 1e-2) && (t<tmax)
            err = oneStep(M,A,VQ[tid]);
            t += 1
            if ParsedARGS["timeseries"]
                h = 1.0 .- min.(VQ[tid].hd, VQ[tid].hu);
                esri = sum([x.sout0 * h[x.id] for x in values(M.Companies)]) / total_volume;
                push!(VR[tid], (i, esri, t));
            end
        end
        if !ParsedARGS["timeseries"]
            h = 1.0 .- min.(VQ[tid].hd, VQ[tid].hu);
            esri = sum([x.sout0 * h[x.id] for x in values(M.Companies)]) / total_volume;
            push!(VR[tid], (i, esri, t));
        end
    end

    df = sort(vcat(VR...), [:index,:t]);
    return df;
end

"""
    saveESRI(M::Market, esri::DataFrame, outfile::String, ParsedARGS)

Saves the computed ESRI results to a CSV file.
- If no custom `psi_mat` was used, it maps indices back to Company names.
- If a custom `psi_mat` was used, it saves Scenario indices.
"""
function saveESRI(M::Market, esri::DataFrame, outfile::String, ParsedARGS)
    if ParsedARGS["psi_mat"] == ""
        companies = getfield.(getindex.(Ref(M.Companies),esri.index),:name)
        dfout = DataFrame(company=companies, esri=esri.esri, t = esri.t);
        CSV.write(outfile, dfout);
    else
        dfout = DataFrame(scenario = esri.index, esri=esri.esri, t=esri.t)
        CSV.write(outfile,dfout)
    end
end