using LinearAlgebra, SparseArrays

abstract type AbstTrjPbm end

"""
General trajectory planning problem class
"""

# General definition and data structure of Auto OPC trajectory generation problem
mutable struct AutoTrjPbm <: AbstTrjPbm

    """ I. Environmental and road paramters """
    # Road info along the guess trajectory


    """ II. Auto Dynamic"""
    # Vehicle model and paramters
    dynmdl::DynMdl

    """ 2.2 Controller constraints """
    # Acceleration(control input) must be slope-limited piecewise continuity
    # 2.2.1 actuators' amplitude limits
    uLowThd::Vector{Float64}     # nu*1 low threshold
    uHighThd::Vector{Float64}    # nu*1 saturation

    # 2.2.3 acceleration(control)'s derivative limits
    uSlopThd::Vector{Float64}    # nu*1

    """ 2.3 State Constraints """
    # States limits generally focuses on speed(first-order derivative)'s amplitude
    # 2-D position states are used for collision-free constraints
    # 2.3.1 vehicle's absolute speed limits
    spdHighThd::Float64
    spdLowThd::Float64      # as soft constraint on elevated road

    """III. Collision-free"""
    # Except 2-D position states, RSS models are also used in collision-free
    """ 3.1 lane and road boundary """
    # 3.1.1 road boundaries without road marking lines

    # 3.1.2 lane boundaries and guidance line

    """ 3.2 Static obstacles """

    """ 3.3 Moving obstacles Game and prediction > drivable boundary """
    
    """ 3.4 RSS risk model """

    """ IV. Boundary Conditions """
    # 4.1 Initial constraint

    # 4.2 Terminal constraint

    # 4.3 Any Intermediate reference point


    """V. Specific traffic problem """
    """ 5.1 Cost function """
    # costRun::Func;      #Running cost
    # costTrmn::Func;     #terminal cost
    
    # costLaneKp::Any;    # lane keep cost in steady ACC

    """ 5.2 Case-specific constraints """
    

    """ VI. Driving preferences and definition  """
    # Definition parameters for trajectory generation of traffic driving
    # stDrvMode::Enum;    #driving style
    # flgEnaEmgyAvd::Bool;   #emergency avoidance without considering traffic rules or some constraints

    """VII. OPC-SCP parameters """
    # Variable Scaling Parameters
                 
    # SCP Algorithm Class
    # cvxSlvrType::ConvexSolverType;      #specify the solver DiscType
    # Convex ScpSubProblem Class
    # cvxSubPbmType::ConvexSubPbmType;

    """ VIII. Comprehensive risk assessment of trajectory """
    # facTrajFeas::double;
    # flgTrajFeas::Bool;
    
    # facDynFeas::double;     # Violation degree of dynamic-related constraints 
    # facRskClls::double;     # each collision risk(negative) and damage level(positive) 
    # rskClls::double;        # all collision risk

    # flgTrajOptm::Bool;
    # timeTraj::double;        # Multi-phase time
    # lossTraj::double;       # Economy: loss and efficiency
    # cmftTraj::double;       # Comfort
    

    """ X. (Shared) Trajectory data-structure"""
    # Initial guess

    # Results
end
