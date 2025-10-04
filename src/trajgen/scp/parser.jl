
""" Transscription and Parser online 
    From  [Trajectory Generation Problem 0&1] and [Discrete Conic Problem 3]
    to SCPPbm and SubPbm Simultaneously
"""
function scp_upd_dyn!(
    subPbm::ScpSubPbm,
    scpPbm::ScpPbm,
    trjPbm::AbstTrjPbm,
    )::Nothing
    
    tstart = Int(time_ns())
    dscrtz!()
    dynDLTV.timeDiscrtz = Int(time_ns() - tstart) / 1e9
    return nothing
end

function scp_upd_!()
end

function scp_upd_uBox!()
end

