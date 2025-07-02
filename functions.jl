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

function buildArrays(M::Market)::Arrays
    # A = Arrays();

    C = M.Companies;
    Edges = M.Edges;

    nrcomp = length(C);

    # set up the weight matrix (fast way)
    fromVect = [l.supplier for l in Edges];
    toVect   = [l.customer for l in Edges];
    weightVect = [l.weight for l in Edges];
    W = sparse(fromVect, toVect, weightVect, nrcomp, nrcomp);

    # W   = spzeros(nrcomp, nrcomp);
    Λd1 = spzeros(nrcomp, nrcomp);
    Λd2 = spzeros(nrcomp, nrcomp);
    Λu  = spzeros(nrcomp, nrcomp);

    sell = sparse(vec(sum(W, dims=2)));
    # buy =  sparse(vec(sum(W, dims=1)));
    for c in values(C)
        c.sout0 = sum([x.weight for x in c.customers]);
        c.sin0 = sum([x.weight for x in c.suppliers]);
    end

    for l in Edges
        from = l.supplier;
        to = l.customer;
        val = W[from, to];
        Λu[to, from] =  val / sell[from];
        # println("$val $from $to $(buy[from]) $(sell[to])")
    end

    ### we build Λd1 and beta
    β = spzeros(nrcomp);
    # id = 2814;
    for id in keys(C)
        company = C[id];
        D = Dict{Int,Vector{Float64}}(); # contains the supplier weights corresponding to the sector
        F = Dict{Int,Vector{Int}}(); # contains the supplier id of the nodes corresponding to the sector
        E = Dict{Int,Int}(); # link type i.e. if the supplier is essential
        nonessential_weight_sum = 0.0;
        for x in company.suppliers # cycle through the suppliers
            suppliernace = x.suppliernace;
            get!(D, suppliernace, Float64[]);
            # append the weight to the corresponding vector addressed by the sector
            push!(D[suppliernace], x.weight);
            get!(F, suppliernace, Int[]);
            # append the id to the corresponding vector addressed by the sector
            push!(F[suppliernace], x.supplier);
            get!(E, suppliernace, x.type);
            nonessential_weight_sum += x.weight; # include ALL links
            # x.type == 0 || (nonessential_weight_sum += x.weight); # exclude useless links
            # x.type == 1 && (nonessential_weight_sum += x.weight); # include non-essential links only
        end
        all = beta = 0.0;
        for sector in keys(D)
            all += sum(D[sector]);
            # if the sector is a necessary one for the node
            if E[sector] == 2
                beta += sum(D[sector]);
                normalize!(D[sector],1);
                for i in eachindex(F[sector])
                    Λd1[F[sector][i], id] = D[sector][i];
                    # println("($(F[sector][i]), $nodeid) = $(D[sector][i])");
                end
            elseif E[sector] == 1
                for i in eachindex(F[sector])
                    Λd2[F[sector][i], id] = D[sector][i] / nonessential_weight_sum;
                end
            end
        end
        if all>0
            β[id] = beta/all;
        end
    end

    # jldsave("data/arrays.jld2"; W, Λu, Λd1, Λd2, β)
    
    return Arrays(Λu, Λd1, Λd2, β);
end

"""
    marketShare(M::Market, Q::DynamicalQuantities)::Nothing

Calculates the market share for each company within its economic sector and updates the marketshare dictionary in the DynamicalQuantities struct `Q`.

This function first computes the total, shock-adjusted sales volume for each sector. It then calculates each company's share of its respective sector's total volume. The function is designed to be numerically robust, preventing division-by-zero errors.

# Arguments
- `M::Market`: The struct containing the entire market structure, including all companies and sectors.
- `Q::DynamicalQuantities`: The struct holding the dynamically evolving quantities, including the `marketshare` dictionary to be updated.

# Returns
- `Nothing`: The function modifies `Q.marketshare` in-place.
"""
function marketShare(M::Market, Q::DynamicalQuantities)::Nothing
    # Alias for convenience to access all companies in the market.
    C = M.Companies;
    
    # Alias for the dictionary that will store the calculated market share for each company.
    # The keys are company IDs and values are their market shares.
    marketshare = Q.marketshare;
    
    # Get the current downstream production level for each company.
    # hd[i] = 1 means full production, hd[i] < 1 means production is reduced by a shock.
    hd = Q.hd;
    
    # Initialize a dictionary to store the total sales volume for each economic sector (by NACE code).
    volumesector = Dict{Int,Float64}();

    # --- Step 1: Calculate the total effective sales volume for each sector. ---
    # Iterate over every company in the market.
    for c in values(C)
        # Calculate the company's shock-adjusted sales (initial sales * current production level).
        # Then, add this value to the total volume of its corresponding sector.
        # `get(volumesector, c.nace, 0.0)` retrieves the current total for the sector, defaulting to 0.0 if the sector hasn't been seen yet.
        volumesector[c.nace] = get(volumesector, c.nace, 0.0) + c.sout0 * hd[c.id];
    end

    # --- Step 2: Calculate the market share for each individual company. ---
    # Iterate over every company again.
    for company in values(C)
        # Get the company's initial, unshocked sales volume.
        sout0 = company.sout0;
        
        # Get the total, shock-adjusted sales volume of the sector this company belongs to.
        vol_sec = get(volumesector, company.nace, 0.0);

        # This check is crucial for numerical stability.
        # It prevents division-by-zero errors if a sector has a total volume of 0.
        if vol_sec > 0
            # If the sector's volume is positive, calculate the market share.
            # The share is the company's initial sales divided by the total sector sales.
            # `min(1.0, ...)` ensures the market share does not exceed 100%.
            marketshare[company.id] = min(1.0, sout0 / vol_sec);
        else
            # If the sector's volume is 0, the company's share in that market is also 0.
            # This handles the `0 / 0` case and prevents `NaN` or `Inf` results.
            marketshare[company.id] = 0.0;
        end
    end

    # The function modifies the marketshare dictionary within Q directly, so there is no explicit return value.
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
        # println(id)
        essentials, non_essentials = downStream(company, A, Q);
        newhd[id] = minimum((essentials, non_essentials, ψ[id]));

        D_u = upStream(company, A, hu);
        newhu[id] = min(D_u, ψ[id]);
    end

    error = max( maximum( abs.(hd .- newhd) ), maximum( abs.(hu .- newhu) ) );

    # garbage collector friendly: copy vectors without changing Q
    Q.hd .= newhd;
    Q.hu .= newhu;

    return error;
end

function ESRI(M::Market, A::Arrays)::Vector{Float64}
    
    nrcomp = length(M.Companies);
    
    @info "Initializing dynamical quantities";
    Q = DynamicalQuantities(nrcomp);

    esri = Vector{Float64}(undef, nrcomp);
    # h    = Vector{Float64}(undef, nrcomp);
    total_volume = sum([x.sout0 for x in values(M.Companies)]);
    u = ones(nrcomp);

    # @showprogress dt=1 desc="Computing..." for t in 1:tmax
    @showprogress for i in sort(collect(keys(M.Companies)))
        Q.psi .= u; Q.psi[i] = 0.0;
        Q.hd .= u; Q.hu .= u;
        err = 1.0;
        while err > 1e-2
            err = oneStep(M,A,Q);
        end
        # println("$i $err");
        h = 1.0 .- min.(Q.hd, Q.hu);
        esri[i] = sum([x.sout0 * h[x.id] for x in values(M.Companies)]) / total_volume;
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