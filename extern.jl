@info "Loading libraries";
using CSV, DataFrames, JLD2, SparseArrays, ProgressMeter;
using LinearAlgebra: normalize!, I;
using ArgParse, SHA;
using PrettyPrint;

"""
    Edge

Represents a directed supply chain link between two firms.

# Fields
- `supplier::Int`: Node ID of the supplier.
- `customer::Int`: Node ID of the customer.
- `suppliernace::Int`: NACE code of the supplier.
- `customernace::Int`: NACE code of the customer.
- `weight::Float64`: The monetary value or volume of the transaction.
- `type::Int`: Categorization of the link for shock propagation:
    - `0`: Do not propagate shocks.
    - `1`: Non-essential input (substitutable).
    - `2`: Essential input (non-substitutable).
"""
struct Edge
    supplier::Int; # neighbor node
    customer::Int;
    suppliernace::Int;
    customernace::Int;
    weight::Float64;
    
    type::Int; #"link's type: 0=do not propagate shocks 1=non-essential 2=essential"
end

"""
    Company

A mutable node representing a single firm in the network.

# Fields
- `name::String`: The identifier or name of the firm (from data source).
- `id::Int`: Internal integer ID assigned during graph construction.
- `nace::Int`: NACE sector code.
- `suppliers::Vector{Edge}`: List of incoming edges (inputs).
- `customers::Vector{Edge}`: List of outgoing edges (sales).
- `sout0::Float64`: Initial total output volume (baseline revenue).
- `sin0::Float64`: Initial total input volume (baseline costs).
"""
mutable struct Company
    name::String; # string with company's name: it's a data's column
    id::Int; # internally assigned ID
    nace::Int; # sector belonging to
    suppliers::Vector{Edge}; # list of links with suppliers
    customers::Vector{Edge}; 
    sout0::Float64; # initial output volume
    sin0::Float64; # initial input volume
end

"""
    Sector

Stores sectoral meta-data, specifically the topology of input essentiality between sectors.

# Fields
- `nace::Int`: The sector's NACE code.
- `essential::Set{Int}`: Set of NACE codes that provide **essential** inputs to this sector.
- `non_essential::Set{Int}`: Set of NACE codes that provide **non-essential** inputs to this sector.
"""
struct Sector
    nace::Int; # sector's code
    # companies::Set{Company}; # set of companies belongind to this sector
    # company_ids::Vector{Int}; # set of companies belongind to this sector
    essential::Set{Int}; # set of essential sectors nace for this sector
    non_essential::Set{Int};
    
    Sector(nace::Int) = new(nace, Set{Int}(), Set{Int}());
end
# Sector(nace::Int) = Sector(nace, Set{Company}(), Vector{Int}(), Set{Int}(), Set{Int}());

"""
    Market

The central container for the economic network topology.

# Fields
- `Companies::Dict{Int, Company}`: Maps internal IDs to `Company` objects.
- `Edges::Vector{Edge}`: A flat list of all connections in the network.
- `Sectors::Dict{Int, Sector}`: Maps NACE codes to `Sector` objects.
- `CompanyID::Dict{String, Int}`: A lookup dictionary mapping original data names to internal integer IDs.
"""
struct Market
    Companies::Dict{Int, Company};
    Edges::Vector{Edge};
    Sectors::Dict{Int, Sector};
    CompanyID::Dict{String, Int}; # maps companies' names to internal IDs

    Market() = new(Dict{Int, Company}(), Edge[], Dict{Int, Sector}(), Dict{String, Int}());
end

"""
    DynamicalQuantities

Holds the state vectors and temporary arrays for the iterative simulation (fixed-point algorithm).

# Fields
- `marketshare::Dict{Int,Float64}`: Dynamic market share of firms within their sector (updates during simulation).
- `hd::Vector{Float64}`: Downstream production level (health) of firms, range [0,1].
- `hu::Vector{Float64}`: Upstream demand level (health) of firms, range [0,1].
- `newhd` / `newhu`: Buffers for storing the next step's state before update.
- `psi::Vector{Float64}`: The exogenous constraint/shock vector for the current scenario.
"""
struct DynamicalQuantities
    marketshare::Dict{Int,Float64};

    hd::Vector{Float64}; # downstream relative production level
    hu::Vector{Float64}; # upstream relative production level
    
    newhd::Vector{Float64}; # downstream relative production level
    newhu::Vector{Float64}; # upstream relative production level
    
    psi::Vector{Float64}; # initial constraint
    
    changed_firms::Vector{Int}; # list of firms that changed hd
    initial_sector_volumes::Dict{Int,Float64}; # 

    function DynamicalQuantities(dim::Int)
        
        marketshare = Dict{Int,Float64}();
        hd = Vector{Float64}(undef, dim);
        hu = Vector{Float64}(undef, dim);
        newhd = Vector{Float64}(undef, dim);
        newhu = Vector{Float64}(undef, dim);
        psi = Vector{Float64}(undef, dim);
    
        changed_firms = Int[];
        sector_volumes = Dict{Int,Float64}();

        return new(marketshare, hd, hu, newhd, newhu, psi, changed_firms, sector_volumes);
    end

end

"""
    Arrays

Stores the sparse matrices required for the linear algebra formulation of the ESRI model.

# Fields
- `lambda_u::SparseMatrixCSC`: The upstream input matrix (Λ_u). Rows are customers, columns are suppliers. Used to pull demand shocks upstream.
- `lambda_d1::SparseMatrixCSC`: The downstream essential input matrix (Λ_d1). Used to push supply shocks downstream for essential goods.
- `lambda_d2::SparseMatrixCSC`: The downstream non-essential input matrix (Λ_d2).
- `beta::SparseVector`: A vector (β) representing the fraction of total inputs that are essential for each firm.
"""
struct Arrays
    lambda_u::SparseMatrixCSC{Float64, Int64};
    lambda_d1::SparseMatrixCSC{Float64, Int64};
    lambda_d2::SparseMatrixCSC{Float64, Int64};
    beta::SparseVector{Float64, Int64};
end