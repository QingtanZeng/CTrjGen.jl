import Pkg
Pkg.activate("./")

## these lines are required for local installations
# Pkg.develop(path="../../SCPToolbox.jl/")
# Pkg.precompile()

using ECOS
using SCPToolbox
using LinearAlgebra
# using Ipopt
# using SCS
# using OSQP

opts = Dict("verbose" => 0)
pbm = ConicProgram(solver = ECOS, solver_options = opts)

q = @new_variable(pbm, "q")

v = @new_variable(pbm, 3, "v")

y = @new_variable(pbm, (2, 5), "y")

p = @new_variable(pbm, "p")
@scale(p, 10, 5);

@scale(y, [2; 4], [3;1])

pbm = ConicProgram(solver = ECOS, solver_options = Dict("verbose" => 0))

x = @new_variable(pbm, 4, "x")
t = @new_variable(pbm, "t")

cstr = @add_constraint(pbm, SOC, "my-soc", (x, t) -> vcat(100*t, x))

l=2; u=9;
A = [Matrix{Float64}(I, 2, 2); -Matrix{Float64}(I, 2, 2);]
b = [-u;-u; l; l]

pbm = ConicProgram();
x = @new_variable(pbm, 2, "x");
cstr = @add_constraint(pbm, NONPOS, "box", x->A*x.+b)

my_pars = Dict("a" => 1, "b" => 2.5)

opts = Dict("verbose" => 0)
pbm = ConicProgram(my_pars; solver = ECOS, solver_options = opts)

x = @new_variable(pbm, 2, "x")
t = @new_variable(pbm, "t")

@add_constraint(pbm, NONPOS, x -> -x*cst["b"])
@add_constraint(pbm, SOC, (x,t) -> vcat(t,x))
@add_constraint(pbm, NONPOS, (x,t) -> x .- [t;t] .+ cst["a"]*cst["b"])

@add_cost(pbm, t -> cst["a"]*t);

exit_status = solve!(pbm)

