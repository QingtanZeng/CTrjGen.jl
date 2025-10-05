
""" Public methods called to SCPPbm and SubPbm"""
function scp_solve!()
    
    while true
        
        # solve ScpSubPbm
        subpbm_solve!()
        
        # save Results from SoluLnrConPgm to ScpSolu
    
        # Calculate defect and Detect Feasibility
        # Update Problem 2&3 from current ScpSolu and TrjPbm
        scp_upd_dyn!()


        # Update SCPPbm and ScpSubPbm buffer from current ScpSolu 
        scp_upd_scppbm!()
        scp_upd_subpbm!()

        # Record to histscp

    end

    # save Results from ScpSolu to TrjPbm


end

function subpbm_solve!()

    # call solver

    # Get intermediate data once for each loop
    solupgm = 
    histpgm = 
    
end

# After Construct, Initialize SCP from guessed trajectory
function scp_init!(subpbm::ScpSubPbm, scppbm::SCPPbm, trjPbm::AbstTrjPbm)::Nothing
    # Initialize SCP-problem2, Sub-problem3, and solver
    
    # SCP-Problem2: transcription
    scp_upd_dyn!()

    # Sub-Problem3: Pre-parse
    scp_init_vc!()          # initial and parse virtual control
    scp_init_dynbox!()      # initial and parse box limits of dynamic system
    scp_init_tr!()          # initial and parse trust region
    scp_init_quacst!()      # initial and parse quadratic cost

end



""" Private methods called from SCPPbm and SubPbm"""
function scp_upd_scppbm!()
    #update reference trajectory from ScpSolu
    scppbm.xref = copy(ScpSolu.xd)
    scppbm.uref = copy(ScpSolu.ud)
    scppbm.pref = copy(ScpSolu.p)

    
end

function scp_upd_subpbm!()
    #update reference trajectory from ScpPbm
    subpbm.xref = scppbm.xref
    
end
    