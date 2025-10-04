using LinearAlgebra, SparseArrays

""" Enumeration definition"""
@enum(CstrtType, ININ, DYN, TMNL, AFF, SOC, TRRG)
@enum(VarsType, STATK, STATKP1, CTRLK, CTRLKP1, PARS, VCP, VCN, LSLK, CSLK, TRRGC)

""" Sub-Data structure of SCP-Problem 2"""
mutable struct ScpParas
    N::Int         # Number of temporal grid nodes
    Nsub::Int      # Number of subinterval integration

    itrScpMax::Int     # Maximum outer loop number
    itrCSlvMax::Int     # Maximum internal loop number

    CSlvType::Int      # Specify default solver type

    epsl_abs::Float64    # Absolute convergence tolerance
    epsl_rel::Float64    # Relative convergence tolerance
    feas_tol::Float64    # Dynamic feasibility tolerance
    q_exit::Float64    # Stopping criterion norm
end
function ScpParas(; N::Int=10, Nsub::Int=10,
    itrScpMax::Int=30, itrCSlvMax::Int=50,
    CSlvType::Int=0,
    epsl_abs::Float64=1e-2, epsl_rel::Float64=3e-3,
    feas_tol::Float64=1e-2)::ScpParas

    q_exit = 1e-2
    prsscp = ScpParas(N, Nsub, itrScpMax, itrCSlvMax,
        CSlvType,
        epsl_abs, epsl_rel,
        feas_tol, q_exit)
    return prsscp
end
mutable struct SCPScaling
    Sx::Matrix{Float64}  # State scaling coefficient matrix
    cx::Vector{Float64}  # State scaling offset vector
    Su::Matrix{Float64}  # Input scaling coefficient matrix
    cu::Vector{Float64}  # Input scaling offset vector
    Sp::Matrix{Float64}  # Parameter scaling coefficient matrix
    cp::Vector{Float64}  # Parameter scaling offset matrix
    iSx::Matrix{Float64} # Inverse of state scaling matrix
    iSu::Matrix{Float64} # Inverse of input scaling matrix
    iSp::Matrix{Float64} # Inverse of parameter scaling coefficient matrix

    function SCPScaling(nx::Int, nu::Int, np::Int)::SCPScaling
        Sx = Float64.(I(nx))
        cx = zeros(Float64, nx)
        Su = Float64.(I(nu))
        cu = zeros(Float64, nu)
        Sp = Float64.(I(np))
        cp = zeros(Float64, np)
        iSx = inv(Sx)
        iSu = inv(Su)
        iSp = inv(Sp)

        scpScl = new(Sx, cx, Su, cu, Sp, cp, iSx, iSu, iSp)
        return scpScl
    end
end
mutable struct ScpSolu          # private data protection
    """ iteritive buffer each sub-problem results """
    # Properties of SCP
    flgFsb::Bool
    flgOpt::Bool
    itrScp::Int
    itrAllCvx::Int
    timescp::Float64   # [s] time spend

    # cost
    cost::Float64

    # Risk/Feasibility assessment
    # dynamics feasibility
    dfctDyn::Vector{Float64}        # N-1 nodes
    flgFsbDynVec::Vector{Bool}      # N-1 nodes
    flgFsbDyn::Bool                 # overall flag

    # Discrete-time Trajectory :x, u, p
    tNodes::Vector{Float64}         # (Shared) by SCPPbm
    xd::Vector{Vector{Float64}}
    ud::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}

    # continuous-time trjPbm
    # xc
    # uc

    function ScpSolu(scpPrs::ScpParas, dynmdl::DynMdl)::ScpSolu
        N = scpPrs.N
        nx, nu, np = dynmdl.nx, dynmdl.nu, dynmdl.np

        # initial Properties
        flgFsb = false
        flgOpt = false
        itrScp = 0
        itrAllCvx = 0
        timescp = 0.0

        # cost
        cost = Inf

        # Risk/Feasibility assessment
        dfctDyn = zeros(Float64, N - 1)
        flgFsbDynVec = zeros(Bool, N - 1)
        flgFsbDyn = false

        # Discrete-time Trajectory :x, u, p
        tNodes = Vector{Float64}()
        xd = [zeros(Float64, nx) for _ in 1:N]
        ud = [zeros(Float64, nu) for _ in 1:N]
        p = [zeros(Float64, np) for _ in 1:N]

        scpsolu = new(flgFsb, flgOpt, itrScp, itrAllCvx, timescp,
            cost,
            dfctDyn, flgFsbDynVec, flgFsbDyn,
            tNodes, xd, ud, p)

        return scpsolu
    end


end
mutable struct ScpHist          # private data protection
end

""" Discrete OPC problem Transcription and Parsed from [Trajectory Generation Problem 0&1]"""
#=  
1. discrete cost with standard linear or quadratic formula from Problem 0&1
2. affine and conic constraits with standard conic classes  from Problem 0&1
3. standard structure of PTR: 
4. configuration, parameters of PTR; Intermediate results and data from ScpSubPbm; History;
5. all SCP methods: Constructer, Parser
=#
mutable struct SCPPbm           # private data protection
    """ SCP parameters """
    scpPrs::ScpParas
    # problem-specific parameters: timeNodes-grid
    tNodes::Vector{Float64}     # nodes between [0,1]

    """ The standard structure for PTR's parsed discrete OPC """
    # Updated reference Trajectory
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Vector{Float64}}
    # Scaling Matrix
    scpScl::SCPScaling

    # 1.1 parsed dynamic system
    dynDLTV::DLTVSys
    idcsDscrtzSys::IdcsDscrtzSys    # the index of system's 1-D P(tau)
    # 1.2 boundaries
    # 1.3 trust-region constraints
    # 1.4 Affine equalities(Equality constraints: Zero cone K=0, Az=b)
    # 1.5 Affine inequalities(inequality constraints, h-Gz<=0)
    # 1.5 SOC constraints(SOC K2, by h - Gz ∈ K)

    # 1.6 linear cost 
    # 1.7 quadratic cost
    # 1.8 virtual control   
    wtr::Vector{Float64}
    wtrp::Float64
    # 1.9 trust region
    wvc::Float64

    """ Properties of standard structure """
    # Abstract index of PTR's parsed discrete OPC
    #n::Int        # number of all linear equalities
    #nl::Int        # number of all linear inequalities
    #nsoc:Int       # number of all soc inequalities

    """ SCP Solution"""
    soluScp::ScpSolu
    histScp::ScpHist

    # internal constructor for custimized data assignment
    function SCPPbm(trjPbm::AbstTrjPbm,
        scpPrs::ScpParas,
        tNodesE::Union{Vector{Float64},Nothing},
        xrefE::Union{Vector{Vector{Float64}},Nothing},
        urefE::Union{Vector{Vector{Float64}},Nothing},
        prefE::Union{Vector{Vector{Float64}},Nothing},
        scpScl::SCPScaling,
        wtr::Vector{Float64},
        wtrp::Float64,
        wvc::Float64,
    )::SCPPbm
        nx, nu, np = trjPbm.dynmdl.nx, trjPbm.dynmdl.nu, trjPbm.dynmdl.np
        N = scpPrs.N

        # problem-specific parameters: timeNodes-grid
        if isnothing(tNodesE)
            tNodes = collect(range(0, 1, length=N))
        else
            tNodes = copy(tNodesE)
        end
        # Updated reference Trajectory(Shared)
        xref = [zeros(Float64, nx) for _ in 1:N]
        if !isnothing(xrefE)
            xref = deepcopy(xrefE)
        end

        uref = [zeros(Float64, nu) for _ in 1:N]
        if !isnothing(urefE)
            uref = deepcopy(urefE)
        end

        pref = [zeros(Float64, np) for _ in 1:N]
        if !isnothing(prefE)
            pref = deepcopy(prefE)
        end

        # standar structure for PTR's parsed discrete OPC
        dynDLTV = DLTVSys(nx, nu, np, N)
        idcsDscrtzSys = IdcsDscrtzSys(trjPbm.dynmdl)

        # penalty weight of virtual control and trust region

        # SCP Solution
        soluScp = ScpSolu(scpPrs, trjPbm.dynmdl)
        histScp = ScpHist()

        # The constructor `new` must be called with all fields
        scppbm = new(scpPrs, tNodes,
            xref, uref, pref,
            scpScl,
            dynDLTV, idcsDscrtzSys,
            wtr, wtrp, wvc,
            soluScp, histScp)
        return scppbm
    end
end

""" Sub-Data structure of Sub-Problem 3 """
struct IdcsLnrConPgm
    # The indices of linear conic program's z,c,A,b,G,h
    # dimenations of variables
    dims_x::Int         # N*nx
    dims_u::Int
    dims_p::Int
    dims_v::Int
    dims_vc::Int
    dims_sml::Int
    dims_chiEta::Int
    dims_chiEtap::Int
    dims_chic::Int
    lgh_z::Int
    # indices of variables z,c
    idx_x::UnitRange{Int}
    idx_u::UnitRange{Int}
    idx_p::UnitRange{Int}
    idx_v::UnitRange{Int}
    idx_vc::UnitRange{Int}
    idx_sml::UnitRange{Int}
    idx_chiEta::UnitRange{Int}
    idx_chiEtap::UnitRange{Int}
    idx_chic::UnitRange{Int}

    # dimenations of constraints
    dims_ININ::Int       # the dimenations of initial boundary
    dims_DYN::Int        # the dimenations of dynamic constraints
    dims_TMNL::Int       # the dimenations of terminal boundary
    num_AFF::Int         # the number of affine equalities
    dims_AFF::Int        # the dimenations of affine equalities
    dims_XBOX::Int
    dims_UBOX::Int
    dims_PBOX::Int
    dims_TRRG::Int
    dims_SOC::Int
    # indices of constraints
end
function IdcsLnrConPgm(scpPbm::SCPPbm, trjPbm::AbstTrjPbm)::IdcsLnrConPgm
    nx, nu, np = trjPbm.nx, trjPbm.nu, trjPbm.np
    N = scpPbm.scpPrs.N

    # Define variables z
    # dynamic system block v1={x;u;p;v+;v-}
    dims_x = nx * N  #  state x
    dims_u = nu * N  #  control u
    dims_p = np  #  parameters p
    dims_v = 2 * nx * (N + 1)  #  virtual control v={v+;v-}
    # all linear (slack) variables block v2={v';s1;...;sml}
    dims_vc =        # vitual control of nonconvex inequalities
        dims_sml = trjPbm       # linear slack for K≤0 to K=0 
    # all trust region and soc slack variables block v3={s1;...;sml}
    dims_chiEta = trjPbm       # trust reion of dynamic
    dims_chiEtap = trjPbm      # trust region of parameters
    dims_chic = trjPbm         # other soc slack variables

    lgh_z = dims_x + dims_u + dims_p
    +2 * dims_v + dims_vc + dims_sml
    dims_chiEta + dims_chiEtap + dims_chic    # z dimenations

    # Define dynamic parts of Ax=b, block c1={X^P_A, X^b_P}
    # 2 boundaries + (N-1) dynamic nodes constraints, idx0 + idxN + idx(1~(N-1))
    # assume the initial and terminal constraints follows same equation with dynamics
    n0 = nx
    nf = nx
    dimsRow_M_P_A = n0 + (N - 1) * nx + nf   # n0+nf+(N-1)nx
    dimsCol_M_P_A = nx * N + nu * N + np
    +2 * nx * (N + 1)   # same as length(v1)
    size_M_P_A = (dimsRow_M_P_A, dimsCol_M_P_A)
    dims_V_P_b = dimsRow_M_P_A

    # create the Sparse A matrix: I, J, V
    dimsSpElem_M_P_A = n0 * nx + n0 * np + 2 * n0     # idx=0, initial constraint
    +(N - 1) * ((nx * nx + nx) + (2 * nx * nu) + (nx * np) + 2 * nx)
    +nf * nx + nf * np + 2 * nf     # idx=N, terminal constraint

    # Define linear equalities parts
    size_Mb_C_A =
        dims_Mb_C_b =
        # Define the linear slack parts


        # Define trust region parts
        # Define other SOC slack parts 

            return Idcs
end
mutable struct ScpPgmParas
end
mutable struct ScpSubSolu
    solupgm::SoluLnrConPgm
end
""" Discrete Conic problem Parsed from [SCP-Problem 2: Discrete OPC problem]"""
#=  
1. standard programming formula {z,c,A,b,G,h} for a specific conic solver

4. configuration, parameters of Conic programmer; Intermediate results and data from programmer; History;
5. all Sub-Problem methods: Constructer, Parser, Update
=#
mutable struct ScpSubPbm        # private data protection

    #problem-specific parameters: timeNodes-grid """
    tNodes::Vector{Float64}     # (Shared) nodes between [0,1]
    # (Shared) Updated reference Trajectory
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Vector{Float64}}

    """ The standard conic problem from PTR's parsed discrete OPC"""
    # The indices of linear conic program's z,c,A,b,G,h
    idcsLnrConPgm::IdcsLnrConPgm
    # The Matrix A, G
    A::Matrix{Float64}
    G::Matrix{Float64}
    # Linear conic problem including Sparse Matrix
    pgmLnrCon::LnrConPgm
    
    # The solution of subPbm
    solusubpbm::ScpSubSolu

    function ScpSubPbm(scpPbm::SCPPbm, trjPbm::AbstTrjPbm)::ScpSubPbm
        #local variables
        scpPrs = scpPbm.scpPrs
        nx, nu, np = trjPbm.dynmdl.nx, trjPbm.dynmdl.nu, trjPbm.dynmdl.np
        N = scpPrs.N

        #problem-specific parameters: timeNodes-grid """
        tNodes::Vector{Float64}     # (Shared) nodes between [0,1]
        # (Shared) Updated reference Trajectory
        xref::Vector{Vector{Float64}}
        uref::Vector{Vector{Float64}}
        pref::Vector{Vector{Float64}}

        # The indices of linear conic program's z,c,A,b,G,h
        idcsLnrConPgm::IdcsLnrConPgm

        dims_z, dims_b,


        # The Matrix A, G
        A = zeros(dims_b, dims_z)
        G = zeros(dims_b, dims_z)
        G::Matrix{Float64}

        # Define linear objective function c    
        I = zeros(Int64, dimsSpElem_M_P_A)
        J = zeros(Int64, dimsSpElem_M_P_A)
        V = zeros(Float64, dimsSpElem_M_P_A)
        # Linear conic problem including Sparse Matrix
        pgmLnrCon::LnrConPgm

        subpbm = ScpSubPbm(scpPbm, trjPbm, IdcsLnrConPgm, A, G, pgmLnrCon)
        return subpbm
    end
end

