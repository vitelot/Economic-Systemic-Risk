include("extern.jl");
include("functions.jl");

function main(file::String)
    @info "Initializing the market according to input file \"$file\"";
    M = initializeMarket(file);

    @info "Building sparse adjacency matrices";
    A = buildArrays(M);

    @info "Initializing dynamical quantities";
    Q = DynamicalQuantities(length(M.Companies));

    @info "Calculating ESRI";
    esri = ESRI(M, A, Q)
end

esri = main("data/complete_edge_list.csv");
# esri = main("data/test_list.csv");

# @profview main("data/complete_edge_list.csv")

# using Plots
# plot(sort(esri, rev=true))



# argmax(esri)
# M = initializeMarket("data/complete_edge_list.csv");
# M.Companies[1432]




# Q.psi[3739] = 0.0; # Q.psi[965] = 0.0

# psi = Q.psi;
# psi[1] = 0.0;
# oneStep(M, A, Q)
# Q.hd'
# Q.hu'

    # C = M.Companies;
    # a = [x.customers for x in values(C)];
    # largestsupplier = sort!(a, by=length, rev=true); a[1][1].supplier
    # af = filter(x-> length(x)==1, a)

    # b = [x.suppliers for x in values(C)];
    # bf = filter(x-> length(x)==1, b)
    # largestcustomer = sort!(b, by=length, rev=true); b[1][1].customer
