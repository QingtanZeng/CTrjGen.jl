include("../mdl/auto.jl")
include("../trjplan/trjpbm.jl")
include("../cvxslv/lnrConPgm.jl")
include("../scp/scppbm.jl")
include("../scp/discrtz.jl")
include("../scp/parser.jl")
include("../scp/scprun.jl")
include("../utils/rk4.jl")

using LinearAlgebra, SparseArrays,ECOS


"""
This script tests the discretization of dynamic system
"""

mutable struct lnrdyn <: AbstTrjPbm
    # 动力学模型
    dynmdl::DynMdl

    # Initial Guess Trajectory
    tf::Float64
    tNodes::Vector{Float64}     # nodes between [0,1]
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Float64}

    function LnrDyn()::lnrdyn
    # x1_dot = x2
    # x2_dot = (-k/m)*x1 - (c/m)*x2 + (1/m)*u + (1/m)*p
    # Continuous Dynamics: x_dot = A*x + B*u + E*p
    # [ x1_dot ] = [   0  ,   1   ] * [ x1 ] + [  0  ] * u + [  0  ] * p
    # [ x2_dot ]   [ -k/m , -c/m  ]   [ x2 ]   [ 1/m ]       [ 1/m ]
        
    m=1;    # [kg]
    k=2;    # [N/m]
    c=0.5;  # [N*s/m]
    f(x,u,p) = begin
        pos, spd = x
        pull = u
        tf = p
        return [spd;
                -(k/m)*pos-(c/m)*spd+pull/m]
    end
    F(x,u,p) = begin
        pos, spd = x
        pull = u
        tf = p
        return tf * [spd;
                    -(k/m)*pos-(c/m)*spd+pull/m]
    end
    A = begin
        pos, spd = x
        pull = u
        tf = p
        return [0        1;
                -k/m  -c/m]
    end
    B = begin
        pos, spd = x
        pull = u
        tf = p
        return [0;
                1/m]
    end
    E = begin
        pos, spd = x
        pull = u
        tf = p
        return f(x,u,p)
    end
    dynmdl = DynMdl(2, 1, 1, f, F, A, B, E)
    tf = 10.0
    
    

        
    end

end

# define a model of dynamic system
# define the problem
lnrdyn = LnrDyn()
nx, nu, np = trjdb.dynmdl.nx, trjdb.dynmdl.nu, trjdb.dynmdl.np

# configure SCP parameters
prsscptpl = (N=11, Nsub=10, itrScpMax=10, itrCSlvMax=50, feas_tol=1e-5)
prsscp = ScpParas(; prsscptpl...)
# Construct SCP problem and its solution
sclscp = SCPScaling(nx, nu, np)
N = prsscp.N
wtr = 1.0
wtrp = 1.0
wvc = 1000.0
scppbm = SCPPbm(trjdb, prsscp,
    nothing, nothing, nothing, nothing,
    sclscp, wtr, wtrp, wvc)

# Construct the sub-problem and its convex solver
subpbm = ScpSubPbm(scppbm, trjdb)

# 2.0 Initialize Guess, SCP-problem, Sub-problem, and solver
#       including scaling and preparse
    # Guess Initial Trajectory: staight line with 100km/h 10.8s
    tf = sqrt(distance^2+(2*widL)^2)/(100*1e3/3600)     # 3s
    trjdb.tNodes = collect(range(0, 1, length=N))

    trjdb.xref = [zeros(Float64, nx) for _ in 1:N]

    trjdb.xref[1] = copy(trjdb.x_0)
    trjdb.xref[end] = copy(trjdb.x_f)
    
    xref_intp = collect(range(trjdb.x_0, trjdb.x_f, N))
    trjdb.xref[2:end-1] = deepcopy(xref_intp[2:end-1])
    theta = -atand(2*widL/distance)/360*2*pi      # 8.53°
    for x in trjdb.xref[2:end-1]
        x[3] = theta
    end

    trjdb.uref = [[27.778; 0;] for _ in 1:N]
    trjdb.uref[1][2] = -1/360*2*pi      # u1 = -1°/s
    trjdb.uref[end][2] = 1/360*2*pi

    trjdb.pref = [tf,]

    #tstart=time_ns();
    scp_init!(subpbm, scppbm, trjdb)
    #plotlpgm(subpbm, "pgm_init")
    plottrj_init = plot_AutoTrjPbm_DubinCar(scppbm.xref)
    #t_init = Int(time_ns() - tstart) / 1e9

# 3.0 iteritive solving loop
    #tstart=time_ns();
    scp_solve!(subpbm, scppbm, trjdb)
    plottrj_result = plot_AutoTrjPbm_DubinCar(scppbm.xref)
    plottrj_real = plot_AutoTrjPbm_DubinCar(scppbm.dynDLTV.xn)
    display(plottrj_init)
    display(plottrj_result)
    display(plottrj_real)
    #t_solve = Int(time_ns() - tstart) / 1e9

    #println("t_init: $t_init")
    #println("t_solve: $t_solve")

    #= online parsing from Problem2 to Problem3
        scp_upd_dynbox!(subpbm, scppbm, trjdb)
        scp_upd_tr!(subpbm, scppbm, trjdb)

        scp_upd_cost!(subpbm, scppbm, trjdb)

        # review the final conic problem's data
        scp_upd_pgm!(subpbm, scppbm, trjdb)
        plotlpgm(subpbm, "pgm")
        

        # solve ScpSubPbm
        subpbm_solve!(subpbm)
        # save Results from ScpSubSolu to ScpSolu
        # Update SCPPbm and ScpSubPbm buffer from current ScpSolu 
        scp_upd_subpbm!(subpbm, scppbm, trjdb)
        #scp_upd_scppbm!()
    
        # Calculate defect and Detect Feasibility
        # Update Problem 2&3 from current ScpSolu and trjdb
        scp_upd_dyn!(subpbm, scppbm, trjdb)
    =#

# 4.0 Record, assessment, Plot
    
#    return nothing;
# end

# main()
