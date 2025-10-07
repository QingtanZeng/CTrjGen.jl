
mutable struct DLTVSys
    # Shared time nodes and reference Trajectory
    tNodes::Vector{Float64}
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Float64}

    # -(xprop_k+1 − Fref(k))= Ak*xk- Inx*xk+1 + Bk−*uk + Bk+*uk+1 + Ek*p
    # -rn = ...

    # Parsed discrete-time state-equation of dynamic system
    xn::Vector{Vector{Float64}}         #x_2 to x_N from xref_k + rk4
    An::Vector{Matrix{Float64}}
    Bkn::Vector{Matrix{Float64}}
    BkP1n::Vector{Matrix{Float64}}
    En::Vector{Matrix{Float64}}
    rn::Vector{Vector{Float64}}

    # Common info
    timeDiscrtz::Float64        # [s] time taken in system's discretization and update

    """
        DLTV(nx, nu, np, nv, N)

    Basic constructor.

    # Arguments
    - `nx`: state dimension.
    - `nu`: input dimension.
    - `np`: parameter dimension.
    - `nv`: virtual control dimension.
    - `N`: the number of discrete time nodes.

    # Returns
    - `dyn`: the dynamics, with empty (undefined) matrices.
    """
    function DLTVSys(nx::Int64, nu::Int64, np::Int64, N::Int64)::DLTVSys        
        # 为向量的每个元素创建一个指定维度的未初始化矩阵
        xn = [Vector{Float64}(undef, nx) for _ in 1:N-1]
        An = [Matrix{Float64}(undef, nx, nx) for _ in 1:N-1]
        Bkn = [Matrix{Float64}(undef, nx, nu) for _ in 1:N-1]
        BkP1n = [Matrix{Float64}(undef, nx, nu) for _ in 1:N-1]
        En = [Matrix{Float64}(undef, nx, np) for _ in 1:N-1]
        rn = [Vector{Float64}(undef, nx) for _ in 1:N-1]
    
        # 初始化其他字段
        tNodes = Vector{Float64}()
        xref = Vector{Vector{Float64}}()
        uref = Vector{Vector{Float64}}()
        pref = Vector{Float64}()
        timeDiscrtz = 0.0
    
        # 使用 new() 创建并返回实例
        return new(tNodes, xref, uref, pref, xn, An, Bkn, BkP1n, En, rn, timeDiscrtz)
    end
end

function DLTVSys_upd!(  node::Int, P::Vector{Float64}, idcs::IdcsDscrtzSys, 
                        dynmdl::DynMdl, dltv::DLTVSys)::Nothing
    nx, nu, np = dynmdl.nx, dynmdl.nu, dynmdl.np
    # reshape back P_k
    x_kP1 = P[idcs.idx_x]
    A_k = reshape(P[idcs.idx_A], (nx, nx))
    B_kN = reshape(P[idcs.idx_Bk], (nx, nu))
    B_kP = reshape(P[idcs.idx_BkP1], (nx, nu))
    E_k = reshape(P[idcs.idx_E], (nx, np))

    # reference point
    xref_k, uref_k, uref_kP1, pref = dltv.xref[node],dltv.uref[node],dltv.uref[node+1], dltv.pref

    # assignment
    copyto!(dltv.xn[node], x_kP1)
    copyto!(dltv.An[node], A_k)
    copyto!(dltv.Bkn[node], B_kN)
    copyto!(dltv.BkP1n[node], B_kP)
    copyto!(dltv.En[node], E_k)
    copyto!(dltv.rn[node], (x_kP1 - (A_k*xref_k+B_kN*uref_k+B_kP*uref_kP1+E_k*pref)))
    return nothing
end
