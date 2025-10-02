#= the discretization of the continuous state equation of dynamic system =#


struct dscrtzSysIdcs
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
    lgh_P::Int;

    function dscrtzSysIdcs(dynmdl::DynMdl)::dscrtzSysIdcs
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
"""
function dscrtz!(ref:: ,pbm::)::Nothing

# initialize the local variable
    P0 = zeros();

end