
""" Transscription and Parser online 
    From  [Trajectory Generation Problem 0&1] and [Discrete Conic Problem 3]
    to SCPPbm and SubPbm Simultaneously
"""
function scp_upd_dyn!( subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm,
    )::Nothing
    
    # online transcription from Problem1 to Problem2
    tstart = Int(tim_ns())

    xref, uref, pref = scppbm.xref, scppbm.uref, scppbm.pref
    dscrtz!(    xref, uref, pref,
                subpbm, scppbm, trjpbm)

    dynDLTV.timeDiscrtz = Int(time_ns() - tstart) / 1e9

    # online parsing from Problem2 to Problem3
    scp_upd_dynAb!(subpbm, scppbm)

    return nothing
end

""" online parsing from Problem2 to Problem3"""
function scp_upd_dynAb!(subpbm::ScpSubPbm, scppbm::ScpPbm,)::Nothing
    N = trjPbm.scpPrs.N
    dltv = scppbm.dynDLTV
    idcs = scppbm.idcsDscrtzSys

    for node = 1:N-1
        idx_zxk = idcs.idx_zx[node]
        idx_zuk = idcs.idx_zu[node]

        idx_bx = idcs.idx_bx[node]


        A()


    end
end

function scp_upd_dynbox!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
end

function scp_upd_tr!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
end

""" Initial parsing from Problem2 to Problem3"""
function scp_init_bc!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
end

function scp_init_vc!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing          # initial and parse virtual control
end

function scp_init_dynbox!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing      # initial and parse box limits of dynamic system
end 

function scp_init_l1cost!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing      # initial and parse L1-norm cost
end

function  scp_init_tr!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing          # initial and parse trust region
end   

