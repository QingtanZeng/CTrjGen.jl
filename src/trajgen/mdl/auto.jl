struct DynMdl
   # Model dimenations
    nx::Int64
    nu::Int64
    np::Int64

    # Dynamics system without normalization
    f::Function
    # Dynamics system with normalization
    F::Function
    
    # Linearized LTI state equation
    # First-order systems of differential equations
    A::Function    # Linearized LTI state equation
    B::Function    # Linearized LTI state equation
    E::Function    # Derivatives of variable parameters
end

function DynMdl_DubinCar()::DynMdl
    # define the model (tau as [0,1])
    # x = [x , y, theta]
    nx = 3;
    # u = [spd, w]
    nu = 2;
    # p = [tf]
    np = 1;
    # dynamic system
    f(x,u,p)=begin
        xCrd, yCrd, theta = x
        spd, w = u
        return [spd*sin(theta); spd*cos(theta); w]
    end
    # time-normalized dynamic system
    F(x, u, p)=begin
        xCrd, yCrd, theta = x
        spd, w = u
        tau, = p
        return tau*[spd*sin(theta); spd*cos(theta); w]
    end
    # Linearized LTI state equation
    A(x,u,p)=begin
        xCrd, yCrd, theta = x
        spd, w = u
        tau, =p
        return tau*[0 0  spd*cos(theta);
                    0 0 -spd*sin(theta);
                    0 0  0             ;]
    end
    B(x, u,p)=begin
        xCrd, yCrd, theta = x
        spd, w = u
        tau, =p
        return tau*[ sin(theta) 0;
                    cos(theta)  0;
                    0           1]
    end
    E(x,u,p)=begin
        return f(x,u,p)
    end

    dynmdl = DynMdl(nx, nu, np, f, F, A, B, E)

    return dynmdl
end

mutable struct DynCstr
    # dynamic control output constraints, derating from lower-level actuators
    # Control amplitude Box Constraints. Each Control MUST have own BoxLimit
    # uL <= u <= uH
    uLowThd::Vector{Float64}     # nu*1 low threshold
    uHighThd::Vector{Float64}    # nu*1 saturation
    # Control Slope-limited Constraints. Each Control MUST have own SlopLimit
    # 0 <= |d/dt*u| <= uSlop
    uSlopThd::Union{Vector{Float64}, Nothing}    # nu*1

    # dynamic states constraints without collision..., only focus on system
    # 0-order states' amplitude Box Constraints, mass,location:
    # Not each state need owm BoxLimit, so formula:
    # xO0L <= (x_Order0=I_O0*x) <= xO0H
    nx_O0::Int64
    I_O0::Matrix{Float64}          # nxO0*nx
    xO0LowThd::Vector{Float64}     # nxO0*1 low threshold
    xO0HighThd::Vector{Float64}    # nxO0*1 high threshold
    # 1-order states' amplitude Box Constraints, velocity:
    # xO1L <= (x_Order1=I_O1*x) <= xO1H
    nx_O1::Int64
    I_O1::Matrix{Float64}          # nxO1*nx
    xO1LowThd::Vector{Float64}     # nxO1*1 low threshold
    xO1HighThd::Vector{Float64}    # nxO1*1 high threshold

    # dynamic parameters constraints
    pLowThd::Vector{Float64}     # np*1 low threshold
    pHighThd::Vector{Float64}    # np*1 high threshold
end

mutable struct DynBdry
    #
end