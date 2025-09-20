"""
General trajectory planning problem class
"""

using LinearAlgebra
using JuMP
using ..utils

export TrjPbm
export TrjPbm_set_dims!,
        TrjPbm_


export DiscretizationType, FOH, IMPULSE

mutable struct TrjPbm
    # Vehicle model and paramters
    nx::Int;
    nu::Int;
    np::Int;
    # Vehicle Dynamics
    f::Func;    # state time derivative
    A::Func;
    B::Func;
    F::Func;

    
    # scaling parameters

    #Others
    mdl::Any    # Problem-specific data structure

end

function TrjPbm_set_dynamics!(
    pbm::TrjPbm,
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


