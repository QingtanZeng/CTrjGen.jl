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
    n::Int;        # Number of variables

    # linear objective function: Min c'z
    c::Vector{Float64};     # Array of size n, cost function weights

    # Equality constraints: Zero cone K=0, Ax=b, from linear constraints with slack variables
    p::Int;                 # Number of equality constraints, = length(b)

    b::Vector{Float64};     # RHS vector of equalities (can be NULL if no equalities are present)
    A::SparseMatrixCSC{Float64, Int64};   # Sparse A matrix data array (column compressed storage) (can be all NULL if no equalities are present)  

    # inequalities constraints {K>=0, SOC K2}, by h - Gx ∈ K
    m::Int;     # Number of inequality constraints, = length(G rows, or h)
    l::Int;     # Dimension of positive orthant
    ncones::Int;    # Number of second order cones
    q::Vector{Int};     # Array of length 'ncones', defines the dimension of each cone
    nex:Int;      # Number of exponential cones, =0, not used in real-time

    h::Vector{Float64};     # Array of size m, RHS vector of cone constraint
    G::SparseMatrixCSC{Float64, Int64};   # Sparse G matrix data array (column compressed storage)

end

function LnrConPgm(

)::LnrConPgm
end

function LnrConPgm_setup(
    lnrConPgm::LnrConPgm,
)::Any
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
    return none;
end