include("../mdl/auto.jl")
include("../trjplan/trjpbm.jl")
include("../cvxslv/lnrConPgm.jl")
include("../scp/dltv.jl")
include("../scp/discrtz.jl")
include("../scp/scppbm.jl")

using LinearAlgebra, SparseArrays


struct AutoTrjPbm_DubinCar <: AbstTrjPbm
    # 动力学模型
    dynmdl::DynMdl
    # 动力学约束
    dyncstr::DynCstr
end
function AutoTrjPbm_DubinCar()::AutoTrjPbm_DubinCar
    dynmdldubin = DynMdl_DubinCar()

    wAbsMax = 5/360*2*pi            # [5ᵒ/s]>[rad/s] 朝向角 角速度
    spdMin = 60*1e3/3600; spdMax = 135*1e3/3600;     # [Km/h]>[m/s] 高速行驶速度区间
    uHighThd = [spdMax, wAbsMax]    # the highest limit of [spd, w]
    uLowThd  = [spdMin, -wAbsMax]         # the lowest limit of [spd, w]

    nx_O0 = 2
    I_O0 = [1 0 0;
            0 0 1]
    xO0HighThd = [2.5*3.75, 5/360*2*pi]     # 横向三车道; 朝向角最大5°
    xO0LowThd  = [0.5*3.75, -5/360*2*pi]    

    pLowThd = [8.0,]    # 300m 135km/h 最短需要8s
    pHighThd = [16.0,]  # 300m 80km/h 最长需要14s

    dyncstr = DynCstr(uLowThd, uHighThd, nothing,
                      nx_O0, I_O0, xO0LowThd, xO0HighThd,   # 0-order state constraints
                      0, zeros(0,dynmdldubin.nu), Float64[], Float64[],  # 1-order state constraints not used
                      pLowThd, pHighThd)
    autotrjpbm = AutoTrjPbm_DubinCar(dynmdldubin, dyncstr)
    return autotrjpbm
end

#function main()::Nothing

# 1.0 configure and Construct all problem and their data-structure

    #Common parameters or precaculated variables

    # define a model of dynamic system
    # define the problem
    trjdb=AutoTrjPbm_DubinCar()

    # configure SCP parameters
    prsscptpl=(N=10, Nsub=10, itrScpMax=30, itrCSlvMax=50)
    prsscp=ScpParas(;prsscptpl...)
    # Construct SCP problem and its solution
    sclscp = SCPScaling(trjdb.dynmdl.nx, trjdb.dynmdl.nu, trjdb.dynmdl.np)
    wtr = ones(Float64, prsscp.N)
    wtrp = 1.0
    wvc = 1.0
    scppbm=SCPPbm(trjdb, prsscp, 
                    nothing, nothing, nothing, nothing,
                    sclscp, wtr, wtrp, wvc)

    # Construct the sub-problem and its convex solver
    subpbm=ScpSubPbm()

# 2.0 Initialize Guess, SCP-problem, Sub-problem, and solver
#       including scaling and preparse
    scp_init!()

# 3.0 iteritive solving loop
    scp_solve!()


# 4.0 Record, assessment, Plot
    
#    return nothing;
# end

# main()
