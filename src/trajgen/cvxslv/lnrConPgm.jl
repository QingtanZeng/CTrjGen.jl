using LinearAlgebra
using SparseArrays
using ECOS


""" definition and data-struct for the interface of linear conic program """
# take ECOS as default solver
#    min c'z 
#    s.t. A*z = b,
#         h - G*z ∈ K 
#    (where K is a composite cone)
mutable struct LnrConPgm
    # variables vector z
    n::Int64;        # Number of variables
    z::Vector{Float64};     # Array of size n, variables

    # linear objective function: Min c'z
    c::Vector{Float64};     # Array of size n, cost function weights

    # Equality constraints: Zero cone K=0, Ax=b, from linear constraints with slack variables
    p::Int64;                 # Number of equality constraints, = length(b)

    b::Vector{Float64};     # RHS vector of equalities (can be NULL if no equalities are present)
    A::SparseMatrixCSC{Float64, Int64};   # Sparse A matrix data array (column compressed storage) (can be all NULL if no equalities are present)  

    # inequalities constraints {K>=0, SOC K2}, by h - Gx ∈ K
    m::Int64;     # Number of inequality constraints, = length(G rows, or h)
    l::Int64;     # Dimension of positive orthant
    ncones::Int64;    # Number of second order cones
    q::Vector{Int};     # Array of length 'ncones', defines the dimension of each cone
    nex::Int64;      # Number of exponential cones, =0, not used in real-time

    h::Vector{Float64};     # Array of size m, RHS vector of cone constraint
    G::SparseMatrixCSC{Float64, Int64};   # Sparse G matrix data array (column compressed storage)

end

mutable struct SoluLnrConPgm
    #  updated states of solver in one iteraion of one subPbm
    z::Vector{Float64}  # primal variables
                        # dual variables
    pcost::Float64
    dcost::Float64
    gap::Float64
    tau::Float64

    # Inner constructor with arguments
    function SoluLnrConPgm()::SoluLnrConPgm
        z = Vector{Float64}()
        pcost = 0.0
        dcost = 0.0
        gap = 0.0
        tau = 0.0
        
       solupgm = new(z, pcost, dcost, gap, tau)
       return solupgm
    end
end

mutable struct HistLnrConPgm
    # history states of solver in one subPbm
    pgm_hist::Vector{SoluLnrConPgm}
end

#=
function LnrConPgm_setup(
                        lnrConPgm::LnrConPgm,
                        )::Nothing
    pwork = ECOS.ECOS_setup(
                lnrConPgm.n,
                lnrConPgm.m,
                lnrConPgm.p,
                lnrConPgm.l,
                lnrConPgm.ncones,
                lnrConPgm.q,
                lnrConPgm.nex,
                lnrConPgm.G.nzval,
                lnrConPgm.G.colptr .-1,
                lnrConPgm.G.rowval .-1,
                lnrConPgm.A.nzval,
                lnrConPgm.A.colptr .-1,
                lnrConPgm.A.rowval .-1,
                lnrConPgm.c,
                lnrConPgm.h,
                lnrConPgm.b, )
    return nothing
end
=#