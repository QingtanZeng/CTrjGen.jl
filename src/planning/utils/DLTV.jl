"""Discerte-time linear time-variant system (DLTV) """

mutable struct DLTV
    # x[:,k+1] = ...
    A::RealTensor     # ...  A[:, :, k]*x[:, k]+ ...
    Bm::RealTensor    # ... +Bm[:, :, k]*u[:, k]+ ...
    Bp::RealTensor    # ... +Bp[:, :, k]*u[:, k+1]+ ...
    F::RealTensor     # ... +F[:, :, k]*p+ ...
    r::RealMatrix     # ... +r[:, k]+ ...
    E::RealTensor     # ... +E[:, :, k]*v
    timing::RealTypes # [s] Time taken to discretize


    """
    DLTV(nx, nu, np, nv, N)
    """
    function DLTV(nx::Int, nu::Int, np,::Int, nv::Int, N::Int)
        A = RealTensor(undef, nx, nx, N-1);
        Bm = RealTensor(undef, nx, nu, N-1);
        Bp = RealTensor(undef, nx, nu, N-1);
        F = RealTensor(undef, nx, np, N-1);
        r = RealMatrix(undef, nx, N-1);
        E = RealTensor(undef, nx, nv, N-1);
        timing = 0.0;   

        dyn = new(A, Bm, Bp, F, r, E, timing);

        return dyn

    end
end