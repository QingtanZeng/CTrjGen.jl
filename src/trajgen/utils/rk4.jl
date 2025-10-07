"""
    4th Runge-Kutta Core Solver Step
"""
function rk4_core_step(
    ddt_F::Function, 
    x::Vector{Float64}, 
    t::Float64, 
    tp::Float64)::Vector{Float64}
    # Calculate the time step
    h = tp-t
    # Calculate the intermediate slopes
    s1=ddt_F(t, x)
    s2=ddt_F(t + h/2, x + h/2 * s1)
    s3=ddt_F(t + h/2, x+ h/2* s2)
    s4=ddt_F(t + h, x + h*s3)

    xp = x + h/6 * (s1+ 2*(s2+s3) + s4)
    return xp    
end

"""
    RK4 Interface
    Nomalized function and duraiton are best
"""
function rk4(
    ddt_F::Function,
    x0::Vector{Float64},
    tgrid::Vector{Float64};
    full::Bool = false,
)::Union{Vector{Float64}, Vector{Vector{Float64}}}

    nx = length(x0)
    N = length(tgrid)
    if full 
        X = [zeros(Float64, nx) for _ in 1:N]
        copyto!(X[1], x0)
    else
        X = zeros(Float64, nx)
        copyto!(X,x0)
    end

    for k =2:N
        if full==true
            copyto!(X[k], rk4_core_step(ddt_F, X[k-1], tgrid[k-1], tgrid[k]))
        else
            copyto!(X, rk4_core_step(ddt_F, X, tgrid[k-1], tgrid[k])) 
        end
    end

    return X;    
end