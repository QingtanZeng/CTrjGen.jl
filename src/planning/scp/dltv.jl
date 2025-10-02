
mutable struct DLTVSys
    # Shared time nodes and reference Trajectory
    tNodes::Vector{Float64}
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Vector{Float64}}

    # Parsed discrete-time state-equation of dynamic system
    An::Vector{Matrix{Float64}}
    Bkn::Vector{Matrix{Float64}}
    BkP1n::Vector{Matrix{Float64}}
    En::Vector{Matrix{Float64}}
    Pn::Vector{Matrix{Float64}}

    # Common info
    timeDiscrtz::Float64        # the time used in system's discretization

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
    function DLTVSys(nx::Int, nu::Int, np::Int, N::Int)::DLTVSys        
        # 为向量的每个元素创建一个指定维度的未初始化矩阵
        An = [Matrix{Float64}(undef, nx, nx) for _ in 1:N]
        Bkn = [Matrix{Float64}(undef, nx, nu) for _ in 1:N]
        BkP1n = [Matrix{Float64}(undef, nx, nu) for _ in 1:N]
        En = [Matrix{Float64}(undef, nx, np) for _ in 1:N]
        Pn = [Matrix{Float64}(undef, nx, 1) for _ in 1:N] # 假设 Pn 是列向量
    
        # 初始化其他字段
        tNodes = Vector{Float64}(undef, N)
        xref = [Vector{Float64}(undef, nx) for _ in 1:N]
        uref = [Vector{Float64}(undef, nu) for _ in 1:N]
        pref = [Vector{Float64}(undef, np) for _ in 1:N]
        timeDiscrtz = 0.0
    
        # 使用 new() 创建并返回实例
        return new(tNodes, xref, uref, pref, An, Bkn, BkP1n, En, Pn, timeDiscrtz)
    end
end