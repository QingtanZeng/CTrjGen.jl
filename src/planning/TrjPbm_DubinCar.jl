include("automdl/AutoMDL.jl")
include("cvxsolver/CVXSOLVER.jl")
include("utils/UTILS.jl")

using LinearAlgebra, SparseArrays


""" Dubin Car as first constructive example for OPC-SCP"""
mutable struct TrjPbm_DubinCar
   """ II. Auto Dynamic"""
    dynmdl::DynMdl
   
   """ 1.1 The model of Dynamic system """

    """ 1.2 Controller constraints """
    # Acceleration(control input) must be slope-limited piecewise continuity
    # 1.2.1 actuators' amplitude limits

    # 1.2.3 acceleration(control)'s derivative limits

    """ 1.3 State Constraints """
    # States limits generally focuses on speed(first-order derivative)'s amplitude
    # 2-D position states are used for collision-free constraints
    # 1.3.1 vehicle's absolute speed limits

    """III. Collision-free"""
    # Except 2-D position states, RSS models are also used in collision-free
    """ 2.1 lane and road boundary """
    # 2.1.1 road boundaries without road marking lines

    # 2.1.2 lane boundaries and guidance line

    """ 2.2 Static obstacles """

    """ 2.3 Moving obstacles Game and prediction > drivable boundary """

    """ IV. Boundary Conditions """
    # 3.1 Initial constraint

    # 3.2 Terminal constraint

    """ X. Trajectory data-structure"""
    # Initial guess

    # Results


    """ """

end

function TrjPbm_DubinCar(;
    nx::Int = 3,
    nu::Int = 2,
    np::Int = 1,
    f::

)::TrjPbm_DubinCar
    
    TrjPbm = TrjPbm_DubinCar()
    return TrjPbm
end






function main()

#Common parameters or precaculated variables
tf = 5;    # [s]

# define a model of dynamic system
dynMdlDubin = DynMdl_DubinCar()

# define the problem(constraints)



end

main()



