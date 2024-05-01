@info "Loading libraries";
using CSV, DataFrames, JLD2, SparseArrays, LinearAlgebra, ProgressMeter;

struct Edge
    supplier::Int; # neighbor node
    customer::Int;
    suppliernace::Int;
    customernace::Int;
    weight::Float64;
    type::Int;
end

mutable struct Company
    name::String; # string with company's name: it's a data's column
    id::Int; # internally assigned ID
    nace::Int; # sector belonging to
    suppliers::Vector{Edge}; # list of links with suppliers
    customers::Vector{Edge}; 
    sout0::Float64; # initial output volume
    sin0::Float64; # initial input volume
end

struct Sector
    nace::Int; # sector's code
    # companies::Set{Company}; # set of companies belongind to this sector
    # company_ids::Vector{Int}; # set of companies belongind to this sector
    essential::Set{Int}; # set of essential sectors nace for this sector
    non_essential::Set{Int};
end
# Sector(nace::Int) = Sector(nace, Set{Company}(), Vector{Int}(), Set{Int}(), Set{Int}());
Sector(nace::Int) = Sector(nace, Set{Int}(), Set{Int}());

struct Market
    Companies::Dict{Int, Company};
    Edges::Vector{Edge};
    Sectors::Dict{Int, Sector};
    CompanyID::Dict{String, Int}; # maps companies' names to internal IDs
end
Market() = Market(Dict{Int, Company}(), Edge[], Dict{Int, Sector}(), Dict{String, Int}());

struct DynamicalQuantities
    marketshare::Dict{Int,Float64};

    hd::Vector{Float64}; # downstream relative production level
    hu::Vector{Float64}; # upstream relative production level
    
    newhd::Vector{Float64}; # downstream relative production level
    newhu::Vector{Float64}; # upstream relative production level
    
    psi::Vector{Float64}; # initial constraint
end
DynamicalQuantities(dim::Int) = DynamicalQuantities(Dict{Int,Float64}(), ones(dim), ones(dim), ones(dim), ones(dim), ones(dim));

struct Arrays
    lambda_u::SparseMatrixCSC{Float64, Int64};
    lambda_d1::SparseMatrixCSC{Float64, Int64};
    lambda_d2::SparseMatrixCSC{Float64, Int64};
    beta::SparseVector{Float64, Int64};
end