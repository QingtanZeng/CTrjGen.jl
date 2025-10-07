using LinearAlgebra
using SparseArrays
using ECOS

mutable struct SoluLnrConPgm
    #  updated states of solver in one iteraion of one subPbm
    z::Vector{Float64}  # primal variables
                        # dual variables
    pcost::Float64
    dcost::Float64
    gap::Float64
    tau::Float64

    exitcode::Int

    # Inner constructor with arguments
    function SoluLnrConPgm(n::Int)::SoluLnrConPgm

        z = zeros(Float64, n)       # the original solution of primal variables
        pcost = 0.0
        dcost = 0.0
        gap = 0.0
        tau = 0.0
        exitcode = 0
        
       solupgm = new(z, pcost, dcost, gap, tau, exitcode)
       return solupgm
    end
end

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
    A::SparseMatrixCSC{Float64, Int64};   # (Shared) Sparse A matrix data array (column compressed storage) (can be all NULL if no equalities are present)  

    # inequalities constraints {K>=0, SOC K2}, by h - Gx ∈ K
    m::Int64;     # Number of inequality constraints, = length(G rows, or h)
    l::Int64;     # Dimension of positive orthant
    ncones::Int64;    # Number of second order cones
    q::Vector{Int};     # Array of length 'ncones', defines the dimension of each cone
    nex::Int64;      # Number of exponential cones, =0, not used in real-time

    h::Vector{Float64};     # Array of size m, RHS vector of cone constraint
    G::SparseMatrixCSC{Float64, Int64};   # (Shared) Sparse G matrix data array (column compressed storage)

    pwork::Ptr{ECOS.pwork}              # the interface structure of ECOS solver

    solupgm::SoluLnrConPgm

    function LnrConPgm( n::Int64,   # number of all variables, =column(A)
                        p::Int64,   # number of all Ax=b constraints, =rows(A)
                        A::SparseMatrixCSC{Float64, Int64},     # complete sparse matrix structure filled with -1

                        m::Int64,   # number of all h-Gx<=0 constraintsa, =rows(G)
                        l::Int64,   
                        ncones::Int64,
                        q::Vector{Int},
                        h::Vector{Float64},
                        G::SparseMatrixCSC{Float64, Int64},
                        )::LnrConPgm

        z = zeros(Float64, n)   # set vector z variables
        c = zeros(Float64, n)   # set vector c cost function weights

        b = zeros(Float64, p)   # set RHS vector of equalities

        nex=0

        pwork = Ptr{ECOS.pwork}()

        solupgm = SoluLnrConPgm(n)

        pgm = new(n, z, c, p, b, A, m, l, ncones, q, nex, h, G, pwork,solupgm)
        return pgm
    end

end

mutable struct HistLnrConPgm
    # history states of solver in one subPbm
    pgm_hist::Vector{SoluLnrConPgm}
end