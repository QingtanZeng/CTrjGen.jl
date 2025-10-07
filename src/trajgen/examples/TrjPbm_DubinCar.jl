include("../mdl/auto.jl")
include("../trjplan/trjpbm.jl")
include("../cvxslv/lnrConPgm.jl")
include("../scp/dltv.jl")
include("../scp/scppbm.jl")
include("../scp/discrtz.jl")
include("../scp/parser.jl")
include("../scp/scprun.jl")

using LinearAlgebra, SparseArrays


struct AutoTrjPbm_DubinCar <: AbstTrjPbm
    # 动力学模型
    dynmdl::DynMdl
    # 动力学约束
    dyncstr::DynCstr
    # 边界值
    A0::Matrix{Float64}
    x_0::Vector{Float64}
    AN::Matrix{Float64}
    x_f::Vector{Float64}

    I_xc::Vector{Float64}

    # Initial Guess Trajectory
    tf::Float64
    tNodes::Vector{Float64}     # nodes between [0,1]
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Float64}

end


function AutoTrjPbm_DubinCar()::AutoTrjPbm_DubinCar
    
    # Env parameters
    widL = 3.75

    # dynamic
    dynmdldubin = DynMdl_DubinCar()
    nx, nu, np = dynmdldubin.nx, dynmdldubin.nu, dynmdldubin.np
    # dynamic constraints
    wAbsMax = 1/360*2*pi            # [ᵒ/s]>[rad/s] 朝向角 角速度
    spdMin = 60*1e3/3600; spdMax = 135*1e3/3600;     # [Km/h]>[m/s] 高速行驶速度区间
    uHighThd = [spdMax, wAbsMax]    # the highest limit of [spd, w]
    uLowThd  = [spdMin, -wAbsMax]         # the lowest limit of [spd, w]

    nx_O0 = 2
    I_O0 = [1 0 0;
            0 0 1]
    xO0HighThd = [2.5*widL, 5/360*2*pi]     # 横向三车道; 朝向角最大5°
    xO0LowThd  = [0.5*widL, -5/360*2*pi]    

    tfref = sqrt(distance^2+(2*widL)^2)/(100*1e3/3600)
    pLowThd = [0.5* tfref,]    # distance m 135km/h 最短需要1.333s
    pHighThd = [2*tfref,]   # distance m 80km/h 最长需要2.333s

    dyncstr = DynCstr(uLowThd, uHighThd, nothing,
                      nx_O0, I_O0, xO0LowThd, xO0HighThd,   # 0-order state constraints
                      0, zeros(0,dynmdldubin.nu), Float64[], Float64[],  # 1-order state constraints not used
                      pLowThd, pHighThd)

    # boundaries Conditions: Inx*x=x0, Inx*x=xf: affine as first, follow DLTVsys
    # 4.1 Initial constraint
    A0 = Float64.(I(nx)); x_0 = [2.5*widL; 0; 0];    
    # 4.2 Terminal constraint
    AN = Float64.(I(nx)); x_f = [0.5*widL; distance; 0]; 

    # Construct the initial trajectory data
    tf = 1
    tNodes = Vector{Float64}()
    xref = Vector{Vector{Float64}}()
    uref = Vector{Vector{Float64}}()
    pref = Vector{Float64}()

    # Cost: keep trajectory straight line
    # cost_x = θ²
    # Mc_X = [0.0 0 0; 0 0 0; 0 0 1]
    # cost_x = wxc*[0 0 1]*x
    wxc = 100.0
    I_xc = [0.0 0 1]
    I_xc = [0.0, 0.0, 1.0]
    autotrjpbm = AutoTrjPbm_DubinCar(dynmdldubin, dyncstr,
                                    A0, x_0, AN, x_f,
                                    I_xc,
                                    tf, tNodes, xref, uref, pref)
    return autotrjpbm
end

#function main()::Nothing
    widL = 3.75
    distance = 100

# 1.0 configure and Construct all problem and their data-structure

    #Common parameters or precaculated variables

    # define a model of dynamic system
    # define the problem
    trjdb=AutoTrjPbm_DubinCar()
    nx, nu, np = trjdb.dynmdl.nx, trjdb.dynmdl.nu, trjdb.dynmdl.np

    # configure SCP parameters
    prsscptpl=(N=10, Nsub=10, itrScpMax=30, itrCSlvMax=50)
    prsscp=ScpParas(;prsscptpl...)
    # Construct SCP problem and its solution
    sclscp = SCPScaling(nx, nu, np)
    N = prsscp.N
    wtr = 1.0
    wtrp = 1.0
    wvc = 1.0
    scppbm=SCPPbm(trjdb, prsscp, 
                    nothing, nothing, nothing, nothing,
                    sclscp, wtr, wtrp, wvc)

    # Construct the sub-problem and its convex solver
    subpbm=ScpSubPbm(scppbm, trjpbm)

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
    trjdb.uref[1][2] = -0.1/360*2*pi      # u1 = -1°/s
    trjdb.uref[end][2] = 0.1/360*2*pi

    trjdb.pref = [tf,]

    scp_init!(subpbm, scppbm, trjdb))

# 3.0 iteritive solving loop
    scp_solve!(subpbm, scppbm, trjdb)


# 4.0 Record, assessment, Plot
    
#    return nothing;
# end

# main()
