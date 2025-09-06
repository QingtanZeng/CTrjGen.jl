import Pkg
Pkg.activate(".")

## these lines are required for local installations
# Pkg.develop(path="../../SCPToolbox.jl/")
# Pkg.precompile()

using SCPToolbox
using JuMP
using ECOS
using Ipopt
# using SCS
# using OSQP
using Printf


# add constraints
lin_ncvx(z, v, pbm) = @add_constraint(pbm, ZERO, 
    (z, v)-> (-2*cst["z̄"][1])*z[1] .+ z[2] .+ cst["z̄"][1]^2 .+ v );
hypln(z, pbm) = @add_constraint(pbm, ZERO, z-> z[2]+0.1*z[1]-0.6);
hlfsp(z, pbm) = @add_constraint(pbm, NONPOS, z-> z[1]-0.85);
tr_cnstr(z, η, pbm) = @add_constraint(pbm, SOC, (z, η)-> vcat(η, z.-cst["z̄"]));
vb_cnstr(v, μ, pbm) = @add_constraint(pbm, L1, (v, μ)-> vcat(μ, v));

# cost function
cost_fun(η, μ, z, pbm) = @add_cost(pbm, (η, μ, z)-> 
    z[2].+cst["wtr"]*η.+cst["wvb"]*μ);

#SCP configuration
params = Dict();
params["solver"] = ECOS;
params["opts"]   = Dict("verbose"=>0,"abstol"=>1e-8);
params["scp-verbose"] = true;

#Initial reference
params["z̄"] = zeros(2);
itr_max = 30;
ϵ_conv = 1e-6;
params["wtr"] = 10;
params["wvb"] = 1e4;

z_hist = [zeros(2) for k = 1:itr_max]
μ_hist = zeros(itr_max)
η_hist = zeros(itr_max)
cost_hist = zeros(itr_max);


function scpSolve_ncvxPbm!(params, z_hist, μ_hist, cost_hist)
    k = 1;
    itr=0;

    while true
        pbm, z, η, μ = construct_subpbm(params);
        flgExit = solve!(pbm);

        z_hist[k], η_hist[k], μ_hist[k] = value(z), value(η)[1], value(μ)[1];
        cost_hist[k] = z_hist[k][2];

        if(k==1)
            z_histLst = zeros(2);
        else
            z_histLst = z_hist[k-1];
        end
        flgConv = max(η_hist[k], μ_hist[k], max(abs.(z_hist[k]-z_histLst)...)) <= ϵ_conv ;

        if params["scp-verbose"]
             @printf("itration %d |%s |z: %5.1e, %5.1e  |TrstRegion: %7.1e |VrtlCtrl: %7.1e |Cost: %.2f |\n",
             k,string(flgExit),z_hist[k][1], z_hist[k][2], η_hist[k],abs(μ_hist[k]),cost_hist[k])
        end

        if (k==itr_max) ||flgConv
            itr = k;
            break
        else
            params["z̄"] = z_hist[k];
            k += 1;
        end
    end

    return itr
end


function construct_subpbm(params)
    
    # Construct the subproblem
    pbm = ConicProgram(params;solver = params["solver"],solver_options = params["opts"])
    
    # Create variablefunction scs
    z = @new_variable(pbm, 2, "z")
    v = @new_variable(pbm, "v")
    η = @new_variable(pbm, "η")
    μ = @new_variable(pbm, "μ")  
        
    # Linearized nonconvex constraint
    lin_ncvx(z, v, pbm)
    
    # Hyperplane constraint
    hypln(z, pbm)
    
    # Halfspace constraint
    hlfsp(z, pbm)
        
    # Trust-region penalty constraint
    tr_cnstr(z, η, pbm)
    
    # Virtual buffer penalty constraint
    vb_cnstr(v, μ, pbm)
    
    # Cost function
    cost_fun(η, μ, z, pbm)
    
    return pbm, z, η, μ
    
end;


itr_conv = scpSolve_ncvxPbm!(params, z_hist, η_hist, μ_hist);