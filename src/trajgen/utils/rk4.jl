"""
    4th Runge-Kutta Core Solver Step
"""
function rk4_core_step(
    Drvt::T.Func, 
    x::RealVector, 
    t::RealValue, 
    tp::RealValue)::RealVector
    # Calculate the time step
   h = tp-t;
   
    # Calculate the intermediate slopes
   s1=Drvt(t, x);
   s2=Drvt(t + h/2, x + h/2 * s1);
   s3=Drvt(t + h/2, x+ h/2* s2);
   s4=Drvt(t + h, x + h*s3);

   xp = x + h/6 * (s1+ 2*(s2+s3) + s4);

   return xp;    
end

"""
    RK4 Interface
    Nomalized function and duraiton are best
"""
function rk4(
    Drvr::T.Func,
    x0::RealVector,
    tspan::RealVector;
    full::Bool = false,
    actions::T.SpecialIntegrationActions = T.SpecialIntegrationActions(undef, 0),
)::Unione{RealVector, RealMatrix}
    X = rk4_generic;
    return X;    
end

"""
    RK4 main function
"""
function rk4_generic(


)

    for k =2:N
        if full
            X[:,k] = rk4_core_step(Drvr, X[:,k-1], tspan[k-1], tspan[k]);
            for act in actions
                X[act[1], k] = act[2](X[act[1], k]);
            end

end