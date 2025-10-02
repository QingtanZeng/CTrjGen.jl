using SparseArrays


""" Enumeration definition"""
@enum(CstrtType, ININ, DYN, TMNL, AFF, SOC, TRRG);
@enum(VarsType, STATK, STATKP1, CTRLK, CTRLKP1, PARS, VCP, VCN, LSLK, CSLK, TRRGC);

""" data structure of SCP and its variants """
mutable struct ScpParas <: SCPParameters
    N::Int;         # Number of temporal grid nodes
    Nsub::Int;      # Number of subinterval integration

    itrOutMax::Int;     # Maximum outer loop number
    itrIntMax::Int;     # Maximum internal loop number

    discMthd::DiscType; # The discretization method

    cvxSlvrType::ConvexSolverType;      # Specify default solver type

    epsl_abs::RealTypes;    # Absolute convergence tolerance
    epsl_rel::RealTypes;    # Relative convergence tolerance
    feas_tol::RealTypes;    # Dynamic feasibility tolerance
end

""" Discrete OPC problem Parsed from Trajectory Generation"""
mutable struct SCPPbm
    """ SCP parameters """
    scpPrs::ScpParas
    """ problem-specific parameters: timeNodes-grid """
    tNodes::Vector{Float64}     # nodes between [0,1], Shared
    
    """ The standard structure for PTR's parsed discrete OPC """
    # Updated reference Trajectory(Shared)
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Vector{Float64}}

    # 1.1 parsed dynamic system
        dynDLTV::DLTV
    # 1.2 boundaries
    # 1.3 trust-region constraints
    # 1.4 Affine equalities(Equality constraints: Zero cone K=0, Az=b)
    # 1.5 Affine inequalities(inequality constraints, h-Gz<=0)
    # 1.5 SOC constraints(SOC K2, by h - Gz ∈ K)

    # 1.6 linear cost 
    # 1.7 quadratic cost
    # 1.8 virtual control   
    wtr::Vector{Float64}
    wtrp::Vector{Float64}
    # 1.9 trust region
    wvc::Vector{Float64}

    """ Properties of standard structure """

    """ SCP Solution"""
    soluScp::ScpSolu
    hstyScp::ScpHsty
end

mutable struct ScpSubPbm
    """ The standard conic problem of solver from discrete OPC""" 
    # Constraints {K=0, K≤0, K2}
    nc::Int;        # number of all linear equalities
    nl::Int;        # number of all linear inequalities
    nsoc:Int;       # number of all soc inequalities

    # Linear conic problem including Sparse Matrix
    pgmLnrCon::LnrConPgm;
    # The indices of linear conic program's z,c,A,b,G,h
    IdcsLnrConPgm::IdcsLnrConPgm
end

""" Construct ScpSubPbm """
function ScpSubPbm(scpPrs::ScpParas, scpPbm::SCPPbm, trjPbm::TrjPbm)::ScpSubPbm
    #local variables
    nx, nu, np = trjPbm.nx, trjPbm.nu, trjPbm.np
    N=scpPrs.N;

    scpPbm.pgmLnrCon.n = dims_z;
    scpPbm.pgmLnrCon.z = zeros(dims_z);

    # Define linear objective function c    

    # Define dynamic parts of Ax=b, block c1={X^P_A, X^b_P}
    # 2 boundaries + (N-1) dynamic nodes constraints, idx0 + idxN + idx(1~(N-1))
    # assume the initial and terminal constraints follows same equation with dynamics
    n0=nx;
    nf=nx;
    dimsRow_M_P_A = n0 + (N - 1)*nx + nf;   # n0+nf+(N-1)nx
    dimsCol_M_P_A = nx * N + nu * N + np 
                    + 2*nx * (N+1);   # same as length(v1)
    size_M_P_A = (dimsRow_M_P_A, dimsCol_M_P_A);      
    dims_V_P_b = dimsRow_M_P_A;

    # create the Sparse A matrix: I, J, V
    dimsSpElem_M_P_A = n0*nx + n0*np + 2*n0     # idx=0, initial constraint
                       +(N-1)*((nx*nx+nx)+(2*nx*nu)+(nx*np)+2*nx)   
                       + nf*nx + nf*np + 2*nf;     # idx=N, terminal constraint
    I=zeros(Int64, dimsSpElem_M_P_A);
    J=zeros(Int64, dimsSpElem_M_P_A);
    V=zeros(Float64, dimsSpElem_M_P_A);

    # idx=0, initial constraint, start{(1,1),(1,),(1,),(1,)}

    # idx=1~(N-1), nodes constraints
    for idx = 0:scpPbm.N

            
    end
    # idx=N, terminal constraint
    
    l


    # Define linear equalities parts
    size_Mb_C_A =
    dims_Mb_C_b = 
    # Define the linear slack parts


    # Define trust region parts
    # Define other SOC slack parts 

    
    pgmLnrCon = LnrConPgm(dims_z,  

end

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
    

    function IdcsLnrConPgm(scpPbm::SCPPbm, trjPbm::TrjPbm)::IdcsLnrConPgm
        nx, nu, np = trjPbm.nx, trjPbm.nu, trjPbm.np
        N=scpPbm.scpPrs.N;

        # Define variables z
        # dynamic system block v1={x;u;p;v+;v-}
        dims_x = nx * N;  #  state x
        dims_u = nu * N;  #  control u
        dims_p = np;  #  parameters p
        dims_v = 2*nx * (N+1);  #  virtual control v={v+;v-}
        # all linear (slack) variables block v2={v';s1;...;sml}
        dims_vc =        # vitual control of nonconvex inequalities
        dims_sml = trjPbm       # linear slack for K≤0 to K=0 
        # all trust region and soc slack variables block v3={s1;...;sml}
        dims_chiEta = trjPbm       # trust reion of dynamic
        dims_chiEtap = trjPbm      # trust region of parameters
        dims_chic = trjPbm         # other soc slack variables

        lgh_z = dims_x + dims_u + dims_p 
                + 2* dims_v + dims_vc + dims_sml
                + dims_chiEta + dims_chiEtap + dims_chic;    # z dimenations


    return Idcs
end

mutable struct ScpSolu
    # Properties of SCP
    flgFsb::Bool
    flgOpt::Bool
    itrScp::Int
    itrAllCvx::Int
    time::Float64   # [s] time spend

    # Properties of terminal CvxSolver
    soluPgm::SoluLnrConPgm

    # cost
    cost::Float64
    
    # Risk/Feasibility assessment

    # Discrete-time Trajectory: tNodes, x, u, p
    tNodes::Vector{Float64}
    xd::Vector{Vector{Float64}}
    ud::Vector{Vector{Float64}}
    p::Vector{Vector{Float64}}

    # continuous-time trjPbm
    xc
    uc
end

mutable struct ScpHsty

end

""" Solve """

function scpsubpbm_solve!()
    
    # solve ScpSubPbm

    # 

end

function scp_upd_dyn!(
    subPbm::ScpSubPbm,
    scpPbm::ScpPbm,
    trjPbm::TrjPbm,
    )::Nothing
    
end
