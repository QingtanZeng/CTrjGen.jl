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

mutable struct LnrDyn <: AbstTrjPbm
    # 动力学模型
    dynmdl::DynMdl

    # Initial Guess Trajectory
    tf::Float64
    tNodes::Vector{Float64}     # nodes between [0,1]
    xref::Vector{Vector{Float64}}
    uref::Vector{Vector{Float64}}
    pref::Vector{Float64}
end
function SpgDmp1Ord()::LnrDyn
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
    A(x,u,p) = begin
        pos, spd = x
        pull = u
        tf = p
        return [0        1;
                -k/m  -c/m]
    end
    B(x,u,p) = begin
        pos, spd = x
        pull = u
        tf = p
        return [0;
                1/m]
    end
    E(x,u,p) = begin
        pos, spd = x
        pull = u
        tf = p
        return f(x,u,p)
    end
    dynmdl = DynMdl(2, 1, 1, f, F, A, B, E)
    # Construct the initial trajectory data
    tf = 1
    tNodes = Vector{Float64}()
    xref = Vector{Vector{Float64}}()
    uref = Vector{Vector{Float64}}()
    pref = Vector{Float64}()
    
    return LnrDyn(dynmdl, tf, tNodes, xref, uref,pref)
end

# define a model of dynamic system
# define the problem
lnrdyn = SpgDmp1Ord()
nx, nu, np = lnrdyn.dynmdl.nx, lnrdyn.dynmdl.nu, lnrdyn.dynmdl.np

# configure SCP parameters
prsscptpl = (N=11, Nsub=10, itrScpMax=10, itrPgmMax=50, feas_tol=1e-5)
prsscp = ScpParas(; prsscptpl...)
# Construct SCP problem and its solution
sclscp = SCPScaling(nx, nu, np)
N = prsscp.N
wtr = 1.0
wtrp = 1.0
wvc = 1000.0
scppbm = SCPPbm(lnrdyn, prsscp,
    nothing, nothing, nothing, nothing,
    sclscp, wtr, wtrp, wvc)

# Construct the sub-problem and its convex solver
subpbm = ScpSubPbm(scppbm, lnrdyn)

# 2.0 Initialize Guess, SCP-problem, Sub-problem, and solver
#       including scaling and preparse
    # Guess Initial Trajectory: staight line with 100km/h 10.8s
    tf = sqrt(distance^2+(2*widL)^2)/(100*1e3/3600)     # 3s
    lnrdyn.tNodes = collect(range(0, 1, length=N))

    lnrdyn.xref = [zeros(Float64, nx) for _ in 1:N]

    lnrdyn.xref[1] = copy(lnrdyn.x_0)
    lnrdyn.xref[end] = copy(lnrdyn.x_f)
    
    xref_intp = collect(range(lnrdyn.x_0, lnrdyn.x_f, N))
    lnrdyn.xref[2:end-1] = deepcopy(xref_intp[2:end-1])
    theta = -atand(2*widL/distance)/360*2*pi      # 8.53°
    for x in lnrdyn.xref[2:end-1]
        x[3] = theta
    end

    lnrdyn.uref = [[27.778; 0;] for _ in 1:N]
    lnrdyn.uref[1][2] = -1/360*2*pi      # u1 = -1°/s
    lnrdyn.uref[end][2] = 1/360*2*pi

    lnrdyn.pref = [tf,]
