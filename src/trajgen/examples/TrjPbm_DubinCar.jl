include("PLANNING.jl")
using LinearAlgebra, SparseArrays


function AutoTrjPbm_DubinCar()::AutoTrjPbm
    dynMdlDubin = DynMdl_DubinCar()
    
    autotrjpbm = AutoTrjPbm(dynMdlDubin)
    return autotrjpbm
end


# function main()::Nothing

# 1.0 configure and Construct all problem and their data-structure

    #Common parameters or precaculated variables
    tf = 5;    # [s]

    # define a model of dynamic system
    # define the problem
    trjdb=TrjPbm_DubinCar()

    # configure SCP parameters
    prsscptpl=(N=10, Nsub=10, itrScpMax=30, itrCSlvMax=50)
    prsscp=ScpParas(;prsscptpl...)
    # Construct SCP problem and its solution
    scppbm=SCPPbm(prsscp, trjdb)

    # Construct the sub-problem and its convex solver
    subpbm=ScpSubPbm()

# 2.0 Initialize Guess, SCP-problem, Sub-problem, and solver
#       including scaling and preparse
    scp_init!()

# 3.0 iteritive solving loop
    scp_solve!()


# 4.0 Record, assessment, Plot
    
#    return nothing;
# end

# main()



