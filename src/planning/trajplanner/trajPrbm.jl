"""
General trajectory planning problem class
"""

using LinearAlgebra
using JuMP
using ..utils

export TrjPrblm
export TrjPrblm_set_dims!,
        TrjPrblm_

export DiscretizationType, FOH, IMPULSE

# General data structure of Auto OPC trajectory generation problem

mutable struct TrjPrbm
    """I. auto dynamic"""
    # Vehicle model and paramters
    nx::Int;
    nu::Int;
    np::Int;

    """ 1.1 Vehicle Dynamics """
    # auto model alternatives: 1) 3+ Dof single-track mdl 2) 8-14 Dof double-track mdl 
    f::Func;    # state time derivative
    A::Func;
    B::Func;
    F::Func;

    """ 1.2 Controller constraints """
    # actuators' limits

    """ 1.3 State Constraints """
    # vehicle states' limit, only considering physical boundary and legal boundary
    # Vehicle speed < 135Km/h
    # Vehicle acceleration limit
    # The acceleration is continuous and the slope is limited.



    """II. collision-free"""
    """ 2.1 road boundary """


    """ 2.2 Static obstacles """


    """2.3 Moving obstacles Game and prediction > drivable boundary"""


    
    # scaling parameters

    #Others
    mdl::Any    # Problem-specific data structure

end

function TrjPrblm_set_dynamics!(
    pbm::TrjPrblm,
    f::Func, A::Func, B::Func, F::Func,
)::Nothing
    pbm.f = (t, k, x, u, p) -> f(t, k, x, u, p, pbm);
    pbm.A = !isnothing(A) ? (t, k, x, u, p) -> A(t, k, x, u, p, pbm) :
                            (t, k, x, u, p) -> zeros(pbm.nx, pbm,nx);
    pbm.B =
        !isnothing(A) ? (t, k, x, u, p) -> B(t, k, x, u, p, pbm) :
        (t, k, x, u, p) -> zeros(pbm.nx, pbm.nu)
    pbm.F =
        !isnothing(F) ? (t, k, x, u, p) -> F(t, k, x, u, p, pbm) :
        (x, u, p) -> zeros(pbm.nx, pbm.nu)
    return nothing
end


