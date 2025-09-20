"""
General trajectory planning problem class
"""

using LinearAlgebra
using JuMP
using ..utils

export TrjPbm
export TrjPbm_set_dims!,
        TrjPbm_

export DiscretizationType, FOH, IMPULSE

# General definition and data structure of Auto OPC trajectory generation problem
mutable struct TrjPbm

    """ I. Environmental and road paramters """
    # Road info along the guess trajectory


    """ II. Auto Dynamic"""
    # Vehicle model and paramters
    nx::Int;
    nu::Int;
    np::Int;

    """ 1.1 Vehicle Dynamics """
    # auto model alternatives: 1) 3+ Dof single-track mdl 2) 8-14 Dof double-track mdl 
    f::Func;    # First-order systems of differential equations
    A::Func;    # Linearized LTI state equation
    B::Func;    # Linearized LTI state equation
    F::Func;    # Derivatives of variable parameters

    """ 1.2 Controller constraints """
    # Acceleration(control input) must be slope-limited piecewise continuity
    # 1.2.1 actuators' amplitude limits

    # 1.2.2 actuators' direction limits

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
    
    """ 2.4 RSS risk model """

    """ IV. Boundary Conditions """
    # 3.1 Initial constraint

    # 3.2 Terminal constraint


    """V. Specific traffic problem """
    mdl::Any;           #Problem-specific data structure
    """ 5.1 Cost function """
    costRun::Func;      #Running cost
    costTrmn::Func;     #terminal cost
    
    costLaneKp::Any;    # lane keep cost in steady ACC

    """ 5.2 Case-specific constraints """

    """ 5.3 All constraints construction"""
    

    """ VI. Driving preferences and definition  """
    # Definition parameters for trajectory generation of traffic driving
    stDrvMode::Enum;    #driving style
    flgEnaEmgyAvd::Bool;   #emergency avoidance without considering traffic rules or some constraints

    """VII. OPC-SCP parameters """
    # Variable Scaling Parameters
    scpScl::SCPScaling;             
    # SCP Algorithm Class
    cvxSlvrType::ConvexSolverType;      #specify the solver DiscType
    # Convex ScpSubProblem Class
    cvxSubPbmType::ConvexSubPbmType;

    """ VIII. Comprehensive risk assessment of trajectory """
    facTrajFeas::double;
    flgTrajFeas::Bool;
    
    facDynFeas::double;     # Violation degree of dynamic-related constraints 
    facRskClls::double;     # each collision risk(negative) and damage level(positive) 
    rskClls::double;        # all collision risk

    flgTrajOptm::Bool;
    timeTraj::double;        # Multi-phase time
    lossTraj::double;       # Economy: loss and efficiency
    cmftTraj::double;       # Comfort
    

    """ X. Trajectory data-structure"""
    # Initial guess

    # Results


    """ """

end

function TrjPbm_set_dynamics!(
    pbm::TrjPbm,
    f::Func, A::Func, B::Func, F::Func,
)::Nothing
    pbm.f = (t, k, x, u, p) -> f(t, k, x, u, p, pbm);
    pbm.A = !isnothing(A) ? (t, k, x, u, p) -> A(t, k, x, u, p, pbm) :
                            (t, k, x, u, p) -> zeros(pbm.nx, pbm,nx);
    pbm.B =
        !isnothing(A) ? (t, k, x, u, p) -> B(t, k, x, u, p, pbm) :
        (t, k, x, u, p) -> zeros(pbm.nx, pbm.nu)
    pbm.F =
        !isnothing(F) ? (t, k, x, u, p) -> F(t, k, x, u, p, pbm) :
        (x, u, p) -> zeros(pbm.nx, pbm.nu)
    return nothing
end

function TrjPbm_set_ctrlLmt!(

)

end


