include("../mdl/auto.jl")
include("../trjplan/trjpbm.jl")
include("../cvxslv/lnrConPgm.jl")
include("../scp/scppbm.jl")
include("../scp/discrtz.jl")
include("../scp/parser.jl")
include("../scp/scprun.jl")
include("../utils/rk4.jl")

using LinearAlgebra, SparseArrays,ECOS, MAT, Plots, Serialization,ProfileView

function plotlpgm(subpbm::ScpSubPbm , name::String)::Nothing
    pgmM = [    subpbm.pgmLnrCon.z'     NaN;
                subpbm.c'               NaN;
                subpbm.A                subpbm.b;  ]
    println("构建的矩阵 M 的维度: ", size(pgmM))
    plotspm(pgmM, name)
    return nothing
end
function plotspm(M_in::Union{Matrix{Float64},Vector{Float64}}, name::String)::Nothing

    if typeof(M_in) == Vector{Float64}
        M = reshape(M_in, length(M_in), 1)
    else
        M=M_in
    end
    matwrite("./trjdb_$name.mat", Dict("$name"=>M) )

    # 用 1 代表 NaN，2 代表 0.0，3 代表非零值。
    display_matrix = zeros(Int8, size(M,1), size(M,2));
    for idx in eachindex(M)
        if isnan(M[idx])==true
            display_matrix[idx] = 1
        elseif iszero(M[idx])==true
            display_matrix[idx] = 2
        elseif isone(M[idx])==true
            display_matrix[idx] = 3
        elseif isone(-M[idx])==true
            display_matrix[idx] = 4
        elseif M[idx] < 0
            display_matrix[idx] = 5
        else 
            display_matrix[idx] = 6
        end
    end
    
    # --- 3. 绘图 ---
    println("正在生成图像...")
    # 1. 定义与分类矩阵对应的 RGBA 颜色
    #    索引 1 -> 略透明的浅黑色 (对应 NaN)
    #    索引 2 -> 非常透明的浅灰色 (对应 0.0)
    #    索引 3 -> 红色 (对应 online parsing)
    #    索引 4 -> 蓝色(对应 1和-1)
    colors = [
        RGBA(0.8, 0.8, 0.8, 0.8),  # 非常透明的浅灰色
        :white,                   # 白色 (完全不透明)
        :blue,                   # 亮蓝色 (完全不透明)
        :purple,                   # 紫色 (完全不透明)
        :red,                    # 红色 (完全不透明)
        :green,                   # 绿色 (完全不透明)
    ]
    gr() # Ensure GR backend is active

    # 自定义坐标轴的标签，让图表更清晰
    widblk = 50
    x_ticks_pos = [1:widblk:size(M, 2); size(M, 2)]
    x_labels = ["$i" for i in x_ticks_pos]
    # 最后一列是 'b'，它的标签需要特殊处理
    x_labels[end]= "$(size(M,2)) (b)"

    y_ticks_pos = [1;2;3:widblk:size(M, 1); size(M, 1)]
    y_labels = ["z'", "c'"]
    append!(y_labels, [ "$i" for i in 3:widblk:size(M,1) ]) # A1,:, A2,:, ...
    push!(y_labels, "$(size(M,1))")

    mapmatrix= heatmap(
        display_matrix,
        c = cgrad(colors, categorical=true), # 使用分类调色板
        colorbar = :none,                    # 分类图例通常不需要颜色条
        aspect_ratio = 1,                    # 保证像素是正方形
        yflip = true,                        # 翻转y轴，使[1,1]在左上角
        #axis = nothing,                      # 隐藏坐标轴刻度
        xticks = (x_ticks_pos, x_labels),
        yticks = (y_ticks_pos, y_labels),
        border = :none,                      # 隐藏边框
        title = "$name $(size(M)): Sparse Matrix Visualization\n(Gray: NaN, White: 0.0, Blue/Purple: ±1, red/green:-+)",
        size = (1000, 1000), # 控制输出图像尺寸
        dpi = 300 # 控制输出图像分辨率
    )
    display(mapmatrix)
    savefig(mapmatrix, "./trjdb_$name.png")

    return nothing
end


mutable struct AutoTrjPbm_DubinCar <: AbstTrjPbm
    # 动力学模型
    dynmdl::DynMdl
    # 动力学约束
    dyncstr::DynCstr
    # 边界值
    A0::Matrix{Float64}
    x_0::Vector{Float64}
    AN::Matrix{Float64}
    x_f::Vector{Float64}

    wxc::Float64
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
    I_xc = [0.0, 0.0, 1.0]
    autotrjpbm = AutoTrjPbm_DubinCar(dynmdldubin, dyncstr,
                                    A0, x_0, AN, x_f,
                                    wxc, I_xc,
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
    prsscptpl=(N=10, Nsub=10, itrScpMax=10, itrCSlvMax=50, feas_tol=1.0)
    prsscp=ScpParas(;prsscptpl...)
    # Construct SCP problem and its solution
    sclscp = SCPScaling(nx, nu, np)
    N = prsscp.N
    wtr = 1.0
    wtrp = 1.0
    wvc = 100.0
    scppbm=SCPPbm(trjdb, prsscp, 
                    nothing, nothing, nothing, nothing,
                    sclscp, wtr, wtrp, wvc)

    # Construct the sub-problem and its convex solver
    subpbm=ScpSubPbm(scppbm, trjdb)

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

    #tstart=time_ns();
    scp_init!(subpbm, scppbm, trjdb)
    plotlpgm(subpbm, "pgm_init")
    #t_init = Int(time_ns() - tstart) / 1e9

# 3.0 iteritive solving loop
    #tstart=time_ns();
    scp_solve!(subpbm, scppbm, trjdb)
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
