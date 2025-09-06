#= the discretization of the continuous state equation of dynamic system =#


struct dscrtzSysIndices
    x::IntRange;
    A::IntRange;
    B::Vector{IntRange};
    F::IntRange;        # Indices for S matrix
    r::IntRange;
    E::IntRange;
    length::Int;

end

"""
    Obtain the discrete-time model(transition matrix) for a continuous-time dynamics system
"""
function dscrtz!(ref:: ,pbm::)::Nothing


    V0 = zeros();

end