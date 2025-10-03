struct DynMdl
   # Model dimenations
    nx::Int
    nu::Int
    np::Int

    # Dynamics system
    f::Function
    
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

    dynmdl = DynMdl(nx, nu, np, f, A, B, E)

    return dynmdl
end

mutable struct DynCstr
    # dynamic control output constraints, derating from lower-level actuators
    # Control amplitude Box Constraints:
    uLowThd::Vector{Float64}     # nu*1 low threshold
    uHighThd::Vector{Float64}    # nu*1 saturation
    # Control Slope-limited Constraints:
    uSlopThd::Vector{Float64}    # nu*1

    # dynamic states constraints without collision..., only focus on system
    # 0-order states' amplitude Box Constraints, mass,location:
    x0LowThd::Vector{Float64}     # nx0*1 low threshold
    x0HighThd::Vector{Float64}    # nx0*1 high threshold
    # 1-order states' amplitude Box Constraints, velocity:
    x1LowThd::Vector{Float64}     # nx1*1 low threshold
    x1HighThd::Vector{Float64}    # nx1*1 high threshold
end

mutable struct DynBdry
    #
end