using LinearAlgebra, SparseArrays

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
    p::Vector{Float64}

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
    pref::Vector{Float64}
    # Scaling Matrix
    scpScl::SCPScaling

    # 1.1 parsed dynamic system
    dynDLTV::DLTVSys
    idcsDscrtzSys::IdcsDscrtzSys    # the index of system's 1-D P(tau)
    # 1.2 boundaries
    A0::Matrix{Float64}     # A0*x1 = x_0
    x_0::Vector{Float64}
    AN::Matrix{Float64}     # AN*xn = x_n
    x_f::Vector{Float64}

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
    function SCPPbm(trjpbm::AbstTrjPbm,
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
    # The dimentations of linear conic program's z,c,A,b,G,h
    dims_z::Int
    dims_b::Int
    num_K0::Int
    Dims_K0::Vector{Int}
    dims_K0::Int
    num_K2::Int
    Dims_K2::Vector{Int}
    dims_K2::Int

    # 1.0 dynamic system block v1={x;u;p;vc;} and A_dyn
    # 1.1 dynamic system
    dims_x::Int         # N*nx,
    dims_u::Int         # N*nu
    dims_p::Int         # N*np
    num_dyn::Int        # (N-1)*nx
    # 1.2 boundaries Conditions
    num_ic::Int        # nic, =size(A0,1)
    num_fc::Int        # nfc, =size(AN,1)
    # 1.3 virtual control variables
    dims_vcdyn::Int     # nx0+(N-1)*nx+nxf, 
                        # usually no vc in initial and terminal boundaries
    dims_svc1::Int      # slack variable, =dims_vcdyn, svc1_i >=|vc_i|_1 
    dims_svc2::Int      # auxiliary variable, =dims_vcdyn, -s2=vc-scv<=0, s2>=0
    dims_svc3::Int      # auxiliary variable, =dims_vcdyn,  s3=vc+scv>=0, s3>=0
    # 1.4 z indices of dynamic system, boundaries Conditions and virtual control
    idx_zx::Function
    idx_zu::Function
    idx_zp::UnitRange{Int}
    idx_zvc::Function
    idx_zvc1::Function
    idx_zvc2::Function
    idx_zvc3::Function
    idx_bx::Function
    idx_bic::UnitRange{Int}
    idx_bfc::UnitRange{Int}
    idx_bvcn::Function
    idx_bvcp::Function

    # 2. all linear (slack) variables block v2={v';s1;...;sml}, Affine equalities and A_Aff
    # 2.1 BOX block for state, control parameters
        dims_sx = 




    idx_chiEta::UnitRange{Int}
    idx_chiEtap::UnitRange{Int}
    idx_chic::UnitRange{Int}



    function IdcsLnrConPgm(scppbm::SCPPbm, trjpbm::AbstTrjPbm)::IdcsLnrConPgm
        nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
        N = scppbm.scpPrs.N
       
        # 1.0 dynamic system block v1={x;u;p;vc;} and A_dyn
        # 1.1 dynamic system
        dims_x = nx * N         #  state x
        dims_u = nu * N         #  control u
        dims_p = np             #  parameters p
        num_dyn = nx * (N-1)    # (N-1) dynamic state equalities
        # 1.2 boundaries Conditions
        num_ic = size(scppbm.A0, 1)
        num_fc = size(scppbm.AN, 1)
        # 1.3 virtual control variables
        dims_vcdyn = (N-1)*nx        # usually no vc in initial and terminal boundaries
        dims_svc1 = dims_vcdyn
        dims_svc2 = dims_vcdyn
        dims_svc3 = dims_vcdyn
        # 1.4 z indices of dynamic system, boundaries Conditions and virtual control
        idx_zx = function(k::Int) 
            return ( ((k-1)*nx+1) : k*nx ) end
        idx_zu = function(k::Int) 
            return ( ((k-1)*nu+1) : k*nu ).+ dims_x  end
        idx_zp = (1:dims_p).+ dims_x .+ dims_u
        idx_zvc = function(k::Int) 
            return ( ((k-1)*nx+1) : k*nx ).+ dims_x .+ dims_u .+ dims_p end
        idx_zsvc1 = (k::Int) -> idx_zvc(k).+dims_vcdyn
        idx_zsvc2 = (k::Int) -> idx_zvc(k).+2*dims_vcdyn
        idx_zsvc3 = (k::Int) -> idx_zvc(k).+3*dims_vcdyn

        idx_bx = function(k::Int) 
            return ( ((k-1)*nx+1) : k*nx ).+ num_ic end
        idx_bic = 1:num_ic
        idx_bfc = (1:num_fc).+ num_ic.+ (N-1)*nx
        idx_bvcn = function(k::Int) 
            return ( ((k-1)*nx+1) : k*nx ).+ dims_vcdyn.+num_ic.+num_fc end
        idx_bvcp = (k::Int)-> idx_bvcn(k).+(N-1)*nx
        end

        # 2. all linear (slack) variables block v2={v';s1;...;sml}, Affine equalities and A_Aff
        # 2.1 BOX block for state, control parameters
        dims_sx = 


        # 3. v3={s1;...;sml}, trust region, SOC and A_SOC
        
        # 3.1 trust region auxiliary block
        dims_chiEta = trjPbm       # trust reion of dynamic
        dims_chiEtap = trjPbm      # trust region of parameters

        # 3.2 quadratic cost auxiliary  block

        # 3.3 SOC auxiliary block
        # dims_chic = trjPbm         # other soc slack variables

        
        # all dimenations of {c,z}, {A, b}, {G, h}
        # Define variable z, c^T*z
        dims_z= dims_x + dims_u + dims_p + dims_v 
                + dims_c
                dims_chiEta + dims_chiEtap + dims_chic    # z dimenations 
        # Define variable b, A*z = b
        num_b
        # Define variable h, h-G*z ∈ K


        idcs = new()
        return idcs
    end
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
    pref::Vector{Float64}

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

    function ScpSubPbm(scppbm::SCPPbm, trjpbm::AbstTrjPbm)::ScpSubPbm
        #local variables
        scpPrs = scppbm.scpPrs
        nx, nu, np = trjPbm.dynmdl.nx, trjPbm.dynmdl.nu, trjPbm.dynmdl.np
        N = scpPrs.N

        #problem-specific parameters: timeNodes-grid """
        tNodes = scppbm.tNodes
        # (Shared) Updated reference Trajectory
        xref=scppbm.xref
        uref=scppbm.uref
        pref=scppbm.pref

        # The standard conic problem from PTR's parsed discrete OPC
        # The indices of linear conic program's z,c,A,b,G,h
        idcs = IdcsLnrConPgm()

        # The Matrix A, G and its sparse form
        A = zeros(idcs.dims_b, idcs.dims_z)
        Asp = spzeros(Float64, Int64, idcs.dims_b, idcs.dims_z)
        G = zeros(idcs.dims_K0+idcs.dims_K2, idcs.dims_z)
        Gsp = spzeros(Float64, Int64, idcs.dims_K0+idcs.dims_K2, idcs.dims_z)

        # Linear conic problem including Sparse Matrix
        # include other matrix z,c,b,h
        pgmLnrCon = LnrConPgm(  idcs.dims_z, 
                                idcs.dims_b,  
                                Asp, 
                                idcs.dims_K0+idcs.dims_K2,
                                idcs.dims_K0,
                                idc.num_K2,
                                idcs.Dims_K2,
                                Gsp
                                )

        # Solution of Sub-Problem
        solusubpbm = ScpSubSolu(idcspgm)

        subpbm = ScpSubPbm(tNodes, xref, uref, pref,
                            idcsLnrConPgm, A, G, pgmLnrCon,
                            solusubpbm)
        return subpbm
    end
end

