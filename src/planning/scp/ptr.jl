using LinearAlgebra
using JuMP
using ..utils

export DiscType, FOH, IMPULSE
@enum(DiscType, FOH, IMPULSE)

""" The PTR algorithm's parameters"""
mutable struct Paras <: SCPParameters

    # Discrete methods and Grid nodes
    N::Int;      # number of temporal grid nodes
    Nsub::Int;      # number of subinterval integration time nodes
    discMthd::DiscType;

    # Optimal Problem's Solver Parameters
    itr_max::Int;   ## Maximum number of iterations
    ϵ_abs::RealTypes;       # Absolute convergence tolerance
    ϵ_rel::RealTypes;       # Relative convergence tolerance
    fsbl_tol::RealTypes;    # Dynamic feasibility tolerance
    solver::Module       # The numerical solver to use for the subproblems
    solver_opts::Dict{String,Any} # Numerical solver options

    # Trust region and virtual control
    wvc::RealTypes;     #Virtual control weights
    wtr::RealTypes;     #Trust region weights
    q_tr::RealTypes      # Trust region norm (possible: 1, 2, 4 (2^2), Inf)
    q_exit::RealTypes    # Stopping criterion norm


end