# ..:: Globals ::..

const T = Types
const RealValue = T.RealTypes
const RealArray = T.RealArray
const RealVector = T.RealVector
const RealMatrix = T.RealMatrix
const Optional = Types.Optional



function rk4_core_step(f::T.Func, s_k::RealVector, t::RealValue, t_kP1::RealValue)
    # Calculate the time step
    dt = t_kP1 - t;

    # Calculate the intermediate slopes
    s1 = func(t, s_k);
    s2 = func(t+dt/2, s_k+dt*s1/2);
    s3 = func(t+dt/2, s_k+dt*s2/2);
    s4 = func(t+dt, s_k+dt*s3);

    s_kPls1 = s_k + dt/6* (s1+2*(s2+s3)+s4);

    return s_kPls1
end
