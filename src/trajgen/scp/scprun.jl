using ECOS, LinearAlgebra,SparseArrays


""" Public methods called to SCPPbm and SubPbm"""
function scp_solve!(subpbm::ScpSubPbm, scppbm::SCPPbm, trjpbm::AbstTrjPbm,)::Int8
    scpPrs = scppbm.scpPrs
    soluscp = scppbm.soluscp
    solupgm = subpbm.pgmLnrCon.solupgm
    
    itrscp = 0
    tstart=time_ns()
    while true
        
        itrscp += 1 
        # online parsing from Problem2 to Problem3
        scp_upd_dynbox!(subpbm, scppbm, trjpbm)
        scp_upd_tr!(subpbm, scppbm, trjpbm)

        scp_upd_cost!(subpbm, scppbm, trjpbm)

        # review the final conic problem's data
        scp_upd_pgm!(subpbm, scppbm, trjpbm)

        # solve ScpSubPbm
        subpbm_solve!(subpbm)
        # save Results from ScpSubSolu to ScpSolu
        # Update SCPPbm and ScpSubPbm buffer from current ScpSolu 
        scp_upd_subpbm!(subpbm, scppbm, trjpbm)
        #scp_upd_scppbm!()
    
        # Calculate defect and Detect Feasibility
        # Update Problem 2&3 from current ScpSolu and TrjPbm
        scp_upd_dyn!(subpbm, scppbm, trjpbm)

        # Record to histscp

        # Stopping criterion
        flgOpt =false
        flgFea = (solupgm.exitcode==ECOS.ECOS_OPTIMAL) && ( soluscp.flgFsbDyn == true)
        if flgFea && flgOpt 
            println("--- END:  $itrscp  Feasible and Optimal ---") 
            println()
            soluscp.codescpexit = 1 
            break
        elseif (itrscp == scpPrs.itrScpMax && flgFea)
            println("--- END: $itrscp  Feasible but Unoptimal ---")
            codescpexit = 2
            break
        elseif(itrscp == scpPrs.itrScpMax && !flgFea)
            println("--- failed ---")
            soluscp.codescpexit = -1
            break
        else
            println("---Iteration $itrscp  Continue: Feasible but Unoptimal ---")
            continue
        end
    end

    soluscp.timescp = Int(time_ns() - tstart) / 1e9
    scppbm.soluscp.itrscp = itrscp

    # save Results from ScpSolu to TrjPbm

    return soluscp.codescpexit

end

function subpbm_solve!(subpbm::ScpSubPbm)::Nothing

    # call solver
    pgm = subpbm.pgmLnrCon
    pwork = pgm.pwork
    solupgm = pgm.solupgm

    solupgm.exitcode = ECOS.ECOS_solve(pwork)

    if solupgm.exitcode != ECOS.ECOS_FATAL
        pwork_loaded = unsafe_load(pwork)
        info = unsafe_load(pwork_loaded.info)

        # --- CRITICAL FIX: Safely copy the solution vector ---
        # The original `unsafe_copyto!` is dangerous and causes memory corruption.
        # The correct way is to create a safe, managed copy of the result.
        copyto!(solupgm.z, unsafe_wrap(Array, pwork_loaded.x, pwork_loaded.n))
        solupgm.gap = info.gap
        solupgm.pcost = info.pcost
    else
        error("Unknown problem in solver .")
    end

    ECOS.ECOS_cleanup(pwork, 0)

    # Get intermediate data once for each loop
    #histpgm = 
    # save Results from SoluLnrConPgm to ScpSubSolu
    
    return nothing
end

# After Construct, Initialize SCP from guessed trajectory
function scp_init!(subpbm::ScpSubPbm, scppbm::SCPPbm, trjpbm::AbstTrjPbm)::Nothing
    # Initialize SCP-problem2, Sub-problem3, and solver
    copyto!( scppbm.xref, trjpbm.xref)
    copyto!( scppbm.uref, trjpbm.uref)
    copyto!( scppbm.pref, trjpbm.pref)
    
    # SCP-Problem2: transcription
    scp_upd_dyn!(subpbm,scppbm,trjpbm)

    # Sub-Problem3: Pre-parse
    scp_init_cost!(subpbm,scppbm,trjpbm)

    scp_init_bc!(subpbm,scppbm,trjpbm)
    scp_init_vc!(subpbm,scppbm,trjpbm)
    scp_init_dynbox!(subpbm,scppbm,trjpbm)      # initial and parse box limits of dynamic system
    scp_init_l1cost!(subpbm,scppbm,trjpbm)      # initial and parse L1-norm cost
    scp_init_tr!(subpbm,scppbm,trjpbm)          # initial and parse trust region
    #scp_init_quacst!()      # initial and parse quadratic cost

    # set sparse A, G,
    scp_init_pgm!(subpbm, scppbm,trjpbm)
    # reset z
    scp_reset_z!(subpbm, scppbm, trjdb)

    return nothing
end

function scp_stopcritera!(subpbm::ScpSubPbm, scppbm::SCPPbm,)::Int

end


""" Private methods called from SCPPbm and SubPbm"""
function scp_upd_subpbm!(subpbm::ScpSubPbm, scppbm::SCPPbm, trjpbm::AbstTrjPbm)::Nothing
    #update reference trajectory from pwork
    idcs=subpbm.idcsLnrConPgm
    pgm=subpbm.pgmLnrCon
    solupgm=subpbm.pgmLnrCon.solupgm
    z = solupgm.z

    # update to scppbm
    for k in 1:scppbm.scpPrs.N
        idx_zxk = idcs.idx_zx(k)
        idx_zuk = idcs.idx_zu(k)
        copyto!(scppbm.xref[k], @view z[idx_zxk])
        copyto!(scppbm.uref[k], @view z[idx_zuk])
    end
    copyto!(scppbm.pref, z[idcs.idcs_z[3]])

    # z to next pgm
    copyto!(pgm.z, z)
    
    return nothing
    
end

function scp_upd_scppbm!(subpbm::ScpSubPbm, scppbm::SCPPbm, trjpbm::AbstTrjPbm)::Nothing
    #update reference trajectory from ScpPbm
    subpbm.xref = scppbm.xref
    
end


