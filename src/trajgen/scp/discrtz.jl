using LinearAlgebra, SparseArrays

#= the discretization of the continuous state equation of dynamic system =#

struct IdcsDscrtzSys
    """ The index of d/dt*P(tau)=
    d/dt*[  x(tau)               [   F;
            P                        A*P;
            P_Bk        =            (B*dltk+A*P_Bk);
            P_BkP1                   (B*dltkP1+A*P_BkP1);
            P_E   ]                  (E+A*P_E)]
    """
    idx_x::UnitRange{Int}
    idx_A::UnitRange{Int}
    idx_Bk::UnitRange{Int}
    idx_BkP1::UnitRange{Int}
    idx_E::UnitRange{Int}
    lgh_P::Int64;

    function IdcsDscrtzSys(dynmdl::DynMdl)::IdcsDscrtzSys
        nx, nu, np = dynmdl.nx, dynmdl.nu, dynmdl.np
        idx_x = (1:nx)
        idx_A = idx_x[end].+ (1:nx*nx)
        idx_Bk = idx_A[end].+ (1:nx*nu)
        idx_BkP1 = idx_Bk[end].+ (1:nx*nu)
        idx_E = idx_BkP1[end].+ (1:nx*np)
        lgh_P = length(idx_x) + length(idx_A) + length(idx_Bk) + length(idx_BkP1) + length(idx_E)

        Idcs = new(idx_x, idx_A, idx_Bk, idx_BkP1, idx_E, lgh_P)

        return Idcs
    end
end

"""
    Obtain the discrete-time model(transition matrix) for a continuous-time dynamics system
    Called at two places: 1) initialization, 2) after solving a sub-problem
    Input: reference trajectory, mdl
    Output: DLTVSys
"""
function dscrtz!(   xref::Vector{Vector{Float64}}, uref::Vector{Vector{Float64}}, 
                    pref::Vector{Float64},
                    subpbm::ScpSubPbm, scppbm::SCPPbm, trjpbm::AbstTrjPbm)::Nothing
    # local variables
    nx, nu, np = trjPbm.dynmdl.nx, trjPbm.dynmdl.nu, trjPbm.dynmdl.np
    N = trjPbm.scpPrs.N
    Nsub = trjPbm.scpPrs.Nsub
    tNodes = scppbm.tNodes
    idcs = scppbm.idcsDscrtzSys
    #Output:
    dynDLTV = scppbm.dynDLTV
    dynDLTV.tNodes = tNodes; dynDLTV.xref = xref; dynDLTV.uref = uref; dynDLTV.pref = pref;

    # initialize P0
    P0 = zeros(idcs.lgh_P)
    P0[idcs.idx_A] = vec(Matrix{Float64}(I,nx,nx)) 

    # dynamic feasibility detection


    # discretization loop
    for node = 1:N
        # initialize P0 's x state as xref_k
        P0[idcs.idx_x] = xref[node]

        # more clear packaging of d/dt*P(tau) = func(tau,x(tau)) used in rk4
        ddt_P = (tau, P) -> derivP_node(tau, P, idcs, 
                                        tNodes[node], tNodes[node+1], 
                                        uref[node], uref[node+1], pref[node], 
                                        trjPbm.dynmdl)
        tSubGrid = collect(range(tNodes[node], tNodes[node+1], length=Nsub))
        P1 = rk4_generic(ddt_P, P0, tSubGrid)

        # update DLTV and pgm.A
        DLTVsys_upd()


        # Calculate defect and feasibility


    end

    return nothing
end

""" General packaging of d/dt*P(tau) for discretization"""
function derivP_node(tau::Float64, P::Vector{Float64}, idcs::IdcsDscrtzSys,
                        tauk::Float64, taukP1::Float64, 
                        uk::Vector{Float64}, ukP1::Vector{Float64}, p::Vector{Float64}, 
                        dynmdl::DynMdl)::Vector{Float64}

    # get Linear interpolation coefficient
    sgmN = (taukP1-tau)/(taukP1-tauk)
    sgmP = (tau-tauk)/(taukP1-tauk)
    # get current control u
    u = sgmN * uk + sgmP * ukP1

    ddt_P = derivP(P, idcs, u, sgmN, sgmP, p, dynmdl)
   
    return ddt_P
end

""" derivative function d/dt*P(tau) = func(P,u,σ⁻,σ⁺,p,mdl)
    Input: P, u, p, dynmdl
    Output: ddt_P
"""
function derivP(P::Vector{Float64}, idcs::IdcsDscrtzSys, 
                u::Vector{Float64}, sgmN::Float64,sgmP::Float64, p::Vector{Float64}, 
                dynmdl::DynMdl)::Vector{Float64}
    nx, nu, np = dynmdl.nx, dynmdl.nu, dynmdl.np

    # reshape back P0
    x = P[idcs.idx_x]
    P_fai = reshape(P[idcs.idx_A], (nx, nx))
    P_faiBN = reshape(P[idcs.idx_Bk], (nx, nu))
    P_faiBP = reshape(P[idcs.idx_BkP1], (nx, nu))
    P_faiE = reshape(P[idcs.idx_E], (nx, np))

    # Linearized LTV matrix
    A_tau = dynmdl.A(x, u, p)
    B_tau = dynmdl.B(x, u, p)
    E_tau = dynmdl.E(x, u, p)

    # d/dt*P(tau)
    F = dynmdl.f(x, u, p)
    P_A = A_tau * P_fai
    P_BN = sgmN * B_tau + A_tau * P_faiBN
    P_BP = sgmP * B_tau + A_tau * P_faiBP
    P_E = E_tau + A_tau * P_faiE

    ddt_P=[
        F;
        vec(P_A);
        vec(P_BN);
        vec(P_BP);
        vec(P_E)
    ]

    return ddt_P
end



