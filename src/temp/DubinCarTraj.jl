import Pkg
Pkg.activate(".")

## these lines are required for local installations
# Pkg.develop(path="../../SCPToolbox.jl/")
# Pkg.precompile()

using SCPToolbox
using PyPlot, Colors, LinearAlgebra

using ECOS


pbm = TrajectoryProblem();
n, m, d = 3, 2, 1
problem_set_dims!(pbm, n, m, d)
t_f = 3
f(t, x, u, p) = begin
    x, y, θ = x;
    v, ω = u;
    return [v*sin(θ); v*cos(θ); ω] * t_f;
end

A(t, x, u, p) = begin
    x, y, θ = x
    v, ω = u
    return [0 0 v*cos(θ);
            0 0 -v*sin(θ);
            0 0 0]*t_f
end

B(t, x, u, p) = begin
    x, y, θ = x
    v, ω = u
    return [sin(θ) 0;
            cos(θ) 0;
            0 1]*t_f
end

F(t, x, u, p) = begin
    return zeros(3, 1)
end;

wrap(func) = (t, k, x, u, p, pbm) -> func(t, x, u, p);
problem_set_dynamics!(pbm, wrap(f), wrap(A), wrap(B), wrap(F));


## car doesn't reverse
x_0 = [1;0;0];
x_f = [-0.9; 1.8; 0];

g_ic(x, p) = x-x_0;
g_tc(x, p) = x-x_f;

H_0(x, p) = I(3);
H_f(x, p) = I(3);

wrap(func) = (x, p, pbm) -> func(x, p);
problem_set_bc!(pbm, :ic, wrap(g_ic), wrap(H_0));
problem_set_bc!(pbm, :tc, wrap(g_tc), wrap(H_f));

c_0 = [-0.1; 1]
r_0 = 1
car_width = 0.1
Δr_0 = car_width/2
E_xy = [1 0 0;0 1 0]

s(t, x, u, p) = [(r_0+Δr_0)^2-(E_xy*x-c_0)'*(E_xy*x-c_0)];
C(t, x, u, p) = reshape(2*E_xy'*(c_0-E_xy*x), 1, 3);

alg = :ptr;
wrap(func) = (t, k, x, u, p, pbm) -> func(t, x, u, p);
problem_set_s!(pbm, alg, wrap(s), wrap(C));

Γ(x, u, p) = u'*u;
wrap(func) = (t, k, x, u, p, pbm) -> func(x, u, p);
problem_set_running_cost!(pbm, alg, wrap(Γ));

state_guess(N) = straightline_interpolate(x_0, x_f, N);
input_guess(N) = straightline_interpolate(zeros(2), zeros(2), N);

problem_set_guess!(pbm, (N, pbm) -> begin
    x = state_guess(N)
    u = input_guess(N)
    p = zeros(1)
    return x, u, p
end);

# PTR Parameters
N, Nsub = 20, 20
iter_max = 50
disc_method = FOH
wvc, wtr = 1e3, 1e0
feas_tol = 5e-3
ε_abs, ε_rel = 1e-5, 1e-3
q_tr = Inf
q_exit = Inf
solver, solver_options = ECOS, Dict("verbose"=>0)

pars = PTR.Parameters(N, Nsub, iter_max, disc_method, wvc, wtr, ε_abs,
                      ε_rel, feas_tol, q_tr, q_exit, solver, solver_options);

if alg == :ptr
    ptr_pbm = PTR.create(pars, pbm)
    sol, history = PTR.solve(ptr_pbm)
end;


include("p3_dubin_plotting.jl");
plot_trajectory()
plt.close() # comment out this line

