" Define basic and composite types used in the projoct. "

using LinearAlgebra
using JuMP


const Optional{T} = Union{T, Nothing};

const IntVector = Vector{Int}
const IntRange = UnitRange{Int}

const RealTypes = Union{Int, Float64};
const RealArray{n} = Union{[Array{T, n} for T in collect(RealTypes)]...};
const RealVector = RealArray{1}
const RealMatrix = RealArray{2}
const RealTensor = RealArray{3}

const Func = Union{Function, Nothing};
