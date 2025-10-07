using LinearAlgebra

""" Transscription and Parser online 
    From  [Trajectory Generation Problem 0&1] and [Discrete Conic Problem 3]
    to SCPPbm and SubPbm Simultaneously
"""
function scp_upd_dyn!(subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm,)::Nothing
    
    # online transcription from Problem1 to Problem2
    tstart = Int(tim_ns())

    xref, uref, pref = scppbm.xref, scppbm.uref, scppbm.pref
    dscrtz!(    xref, uref, pref,
                subpbm, scppbm, trjpbm)

    dynDLTV.timeDiscrtz = Int(time_ns() - tstart) / 1e9

    # online parsing from Problem2 to Problem3
    scp_upd_dynAb!(subpbm, scppbm, trjpbm)

    return nothing
end

""" online parsing from Problem2 to Problem3 dealing with c, A, b and z"""
function scp_upd_dynAb!(subpbm::ScpSubPbm, scppbm::ScpPbm,trjpbm::AbstTrjPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    dltv = scppbm.dynDLTV
    A, b = subpbm.A, subpbm.b

    # -(xprop_k+1 − Fref(k))= Ak*xk- Inx*xk+1 + Bk−*uk + Bk+*uk+1 + Ek*p

    for node = 1:N-1
        idx_zxk = idcs.idx_zx(node)
        idx_zxkP1 = idcs.idx_zx(node+1)
        idx_zuk = idcs.idx_zu(node)
        idx_zukP1 = idcs.idx_zu(node+1)
        idx_zp = idcs.idx_zp

        idx_bxk = idcs.idx_bx(node)

        # A's dynamic part
        A[idx_bxk, idx_zxk] = dltv.An[node]
        A[idx_bxk, idx_zxkP1] = -1*Float64.(I(nx))
        A[idx_bxk, idx_zuk] = dltv.Bkn[node]
        A[idx_bxk, idx_zukP1] = dltv.BkP1n[node]
        A[idx_bxk, idx_zp] = dltv.En[node]
        # b's dynamic part
        b[idx_bxk] = -1*dltv.rn[node]
    end
end

function scp_upd_dynbox!(subpbm::ScpSubPbm, scppbm::ScpPbm,trjpbm::AbstTrjPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    A, b = subpbm.A, subpbm.b

    #[ Ixl  IsMax   0    * [x; sMax; sMin] =[Max; Min]
    #  Ixl  0    -IsMin] 
    b[idcs.idcs_bxmax] = Float64.(kron(ones(N), scppbm.xHighThd))
    b[idcs.idcs_bxmin] = Float64.(kron(ones(N), scppbm.xLowThd))

    b[idcs.idcs_bumax] = Float64.(kron(ones(N), scppbm.uHighThd))
    b[idcs.idcs_bumin] = Float64.(kron(ones(N), scppbm.uLowThd))

    b[idcs.idcs_bpmax] = Float64.(scppbm.pHighThd)
    b[idcs.idcs_bpmin] = Float64.(scppbm.pLowThd)
end

function scp_upd_tr!(subpbm::ScpSubPbm, scppbm::ScpPbm,trjpbm::AbstTrjPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    idcs = subpbm.idcsLnrConPgm
    b = subpbm.b
    xref, uref, pref = scppbm.xref, scppbm.uref, scppbm.pref
    
    b[idcs.idx_bptrn] = pref
    b[idcs.idx_bptrp] = pref

    b[idcs.idx_bxtr] = [x for xk in xref for x in xk]

    b[idcs.idx_butr] = [u for uk in uref for u in uk]
end

function scp_upd_cost!(subpbm::ScpSubPbm, scppbm::ScpPbm,trjpbm::AbstTrjPbm):Nothing
    # update the original weights from Problem1&2
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    pgm = subpbm.pgm
    c = subpbm.c
    wxc, wvc, wtrp, wtr = scppbm.wxc, scppbm.wvc, scppbm.wtrp, scppbm.wtr

    # dynamic states, control, parameter, virtual control
    # auxiliary variables of virtual control svc1, svc2, svc3
    c[idc.idcs_z[5]] = wvc * ones(idcs.dims_vcdyn)       # wvc*sum(svc1) = wvc* 1_(N-1)*nx^T * svc1
    # auxiliary variables of box limits
    # auxiliary variables of l1-norm cost, sxc1, sxc2, sxc3
    c[idc.idcs_z[14]] = wxc * ones(idcs.dims_sxc1)
    # trust region's weights: update from defect
    c[idc.idcs_z[17]] = wtrp*[1; zeros(idcs.dims_p)]
    c[idc.idcs_z[18]] = wtr*kron(ones(N), [1; zeros(nx)])
    c[idc.idcs_z[19]] = wtr*kron(ones(N), [1; zeros(nu)])

    
end

function scp_reset_z!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    A, b = subpbm.A, subpbm.b
    
    z = subpbm.pgmLnrCon.z
    dfctDyn = scppbm.soluscp.dfctDyn

    # x, u, p : reference trajectory
    z[idcs.idcs_z[1]] = [x for xk in scppbm.xref for x in xk]
    z[idcs.idcs_z[2]] = [u for uk in scppbm.uref for u in uk]
    z[idcs.idcs_z[3]] = scppbm.pref

    # vc=defect, reset virtual control at beginning, but l1-norm penalty, same as
    # 
    z[idcs.idcs_z[4]] = [vc for vck in dfctDyn for vc in vck]
    z[idcs.idcs_z[5]] = abs.(z[idcs.idcs_z[4]])     # svc1>=|vc|
    z[idcs.idcs_z[6]] = -1*z[idcs.idcs_z[4]] + z[idcs.idcs_z[5]]          # v'-s1+s2=0
    z[idcs.idcs_z[7]] = z[idcs.idcs_z[4]] + z[idcs.idcs_z[5]]          # v'+s1-s3=0

    # auxiliary variables of box limits
    # I_xl*x + s_xmax = x_high => s_xmax = x_high - I_xl*x
    # I_xl*x - s_xmin = x_low  => s_xmin = I_xl*x - x_low
    # These slack variables must be non-negative (s_xmax >= 0, s_xmin >= 0).
    # We initialize them based on the reference trajectory. If the reference violates
    # the bounds, the initial slack variable will be negative, which is acceptable
    xl = kron(I(N), scppbm.I_xl) * z[idcs.idcs_zx(1)]
    z[idcs.idcs_z[8]] = b[idcs.idcs_bxmax] - xl
    z[idcs.idcs_z[9]] = -1*b[idcs.idcs_bxmin] + xl
    # For u
    z[idcs.idcs_z[10]] = b[idcs.idcs_bumax] - kron(I(N), scppbm.I_ul)*z[idcs.idcs_zx(2)]
    z[idcs.idcs_z[11]] =  -1*b[idcs.idcs_bumin] + kron(I(N), scppbm.I_ul)*z[idcs.idcs_zx(2)]
    # For p
    z[idcs.idcs_z[12]] = b[idcs.idcs_bpmax] - scppbm.I_sp*z[idcs.idcs_zx(3)]
    z[idcs.idcs_z[13]] = -1*b[idcs.idcs_bpmin] + scppbm.I_sp*z[idcs.idcs_zx(3)]

    # auxiliary variables of l1-norm cost
    xc = (kron(I(N), scppbm.I_xc)) * z[idcs.idcs_zx[1]]
    z[idcs.idcs_z[14]] = abs.(xc)  #sxc1>=|theta|
    #theta + sxc1 - sxc2 = 0
    z[idcs.idcs_z[15]] = xc + z[idcs.idcs_z[14]]
    #theta - sxc1 + sxc3 = 0
    z[idcs.idcs_z[16]] = -1*xc + z[idcs.idcs_z[14]]

    # auxiliary variables of trust region
    z[idcs.idcs_z[17]] = zeros(idcs.dims_chiptr)    # at beginning p-pref=0
    z[idcs.idcs_z[18]] = zeros(idcs.dims_chixtr)    # same as x, u
    z[idcs.idcs_z[19]] = zeros(idcs.dims_chiutr)

end

function scp_upd_pgm!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
    # reset the defined conic problem and intermediate variables
    idcs = subpbm.idcsLnrConPgm
    pgm = subpbm.pgm
    Achk = subpbm.Achk
    A = subpbm.A
    Asp = subpbm.Asp
    Gsp = subpbm.Gsp

    # reset z
    scp_reset_z!(subpbm, scppbm)

    # check the {z,c,A,b,G,h}
    pgm.c = copy(subpbm.c)
    pgm.b = copy(subpbm.b)
    # update Asp
    for i  in 1:length(Asp.nzval)
        Asp[i] = A[I[i], J[i]]
    end
    pgm.A = Asp
    pgm.G = Gsp

    # reset solver
    LnrConPgm_upd!(subpbm)
end


""" Initial parsing from Problem2 to Problem3"""

function scp_init_cost!(subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    pgm = subpbm.pgm
    c = subpbm.c
    wxc, wvc, wtrp, wtr = scppbm.wxc, scppbm.wvc, scppbm.wtrp, scppbm.wtr

    # dynamic states, control, parameter, virtual control
    c[idc.idcs_z[1]] = kron(ones(N),zeros(nx))
    c[idc.idcs_z[2]] = kron(ones(N),zeros(nu))
    c[idc.idcs_z[3]] = zeros(np)
    c[idc.idcs_z[4]] = kron(ones(N-1), zeros(nx))
    # auxiliary variables of virtual control svc1, svc2, svc3
    c[idc.idcs_z[5]] = wvc * ones(idcs.dims_vcdyn)       # wvc*sum(svc1) = wvc* 1_(N-1)*nx^T * svc1
    c[idc.idcs_z[6]] = zeros(idcs.dims_vcdyn)
    c[idc.idcs_z[7]] = zeros(idcs.dims_vcdyn)
    # auxiliary variables of box limits
    c[idc.idcs_z[8]] = zeros(idcs.dims_sxmax)
    c[idc.idcs_z[9]] = zeros(idcs.dims_sxmin)
    c[idc.idcs_z[10]] = zeros(idcs.dims_sumax)
    c[idc.idcs_z[11]] = zeros(idcs.dims_sumin)
    c[idc.idcs_z[12]] = zeros(idcs.dims_spmax)
    c[idc.idcs_z[13]] = zeros(idcs.dims_spmin)
    # auxiliary variables of l1-norm cost, sxc1, sxc2, sxc3
    c[idc.idcs_z[14]] =  wxc * ones(idcs.dims_sxc1)
    c[idc.idcs_z[15]] = zeros(idcs.dims_sxc2)
    c[idc.idcs_z[16]] = zeros(idcs.dims_sxc3)
    # trust region's weights from setting and defect
    c[idc.idcs_z[17]] = wtrp*[1; zeros(idcs.dims_chiptr-1)]
    c[idc.idcs_z[18]] = wtr*kron(ones(N), [1; zeros(nx)])
    c[idc.idcs_z[19]] = wtr*kron(ones(N), [1; zeros(nu)])

end

function scp_init_bc!(subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    idcs = subpbm.idcsLnrConPgm
    A = subpbm.A
    b = subpbm.b

    A[idcs.idx_bic, idcs.idx_zx(1)] = trjpbm.A0
    b[idcs.idx_bic] = trjpbm.x_0
    A[idcs.idx_bfc, idcs.idx_zx(1)] = trjpbm.AN
    b[idcs.idx_bfc] = trjpbm.x_N

end

function scp_init_vc!(subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm)::Nothing          # initial and parse virtual control
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    A, b = subpbm.A, subpbm.b
    
    # set the (idx_bx, idx_zvc)
    A[idcs.idcs_bic, idcs.idc_z[4]] = zeros(num_ic, dims_vcdyn)
    A[idcs.idcs_b[2], idcs.idc_z[4]] = Float64.(kron(I(N-1), I(nx)))
    A[idcs.idcs_bfc, idcs.idc_z[4]] = zeros(num_fc, dims_vcdyn)

    # set the (idx_bvcn, idx_zsvc1_zsvc2)   v'-s1+s2=0
    A[idcs.idcs_bvcn, idcs.idc_z[4]] = Float64.(kron(I(N-1), I(nx)))    # set the (idx_bvcn_bvcp, idx_zvc)
    A[idcs.idcs_bvcn, idcs.idc_z[5]] = -1*Float64.(kron(I(N-1), I(nx)))
    A[idcs.idcs_bvcn, idcs.idc_z[6]] = Float64.(kron(I(N-1), I(nx)))
    # set the (idx_bvcn, idx_zsvc1_zsvc2)   v'+s1-s3=0
    A[idcs.idcs_bvcp, idcs.idc_z[4]] = Float64.(kron(I(N-1), I(nx)))    # set the (idx_bvcn_bvcp, idx_zvc)
    A[idcs.idcs_bvcp, idcs.idc_z[5]] = Float64.(kron(I(N-1), I(nx)))
    A[idcs.idcs_bvcp, idcs.idc_z[7]] = -1*Float64.(kron(I(N-1), I(nx)))

    # set the b[idx_bvcn_bvcp]
    b[idcs.idcs_bvcn] = zeros(idcs.dims_vcdyn)
    b[idcs.idcs_bvcp] = zeros(idcs.dims_vcdyn)
end

function scp_init_dynbox!(subpbm::ScpSubPbm, scppbm::ScpPbm, trjpbm::AbstTrjPbm)::Nothing      # initial and parse box limits of dynamic system
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    A, b = subpbm.A, subpbm.b

    #[ Ixl  IsMax   0    * [x; sMax; sMin] =[Max; Min]
    #  Ixl  0    -IsMin] 
    A[idcs.idcs_bxmax, idcs.idcs_z[1]] = Float64.(kron(I(N), scppbm.I_xl))
    A[idcs.idcs_bxmax, idcs.idcs_z[8]] = Float64.(kron(I(N), I(idcs.n_sx)))
    b[idcs.idcs_bxmax] = Float64.(kron(ones(N), scppbm.xHighThd))
    A[idcs.idcs_bxmin, idcs.idcs_z[1]] = Float64.(kron(I(N), scppbm.I_xl))
    A[idcs.idcs_bxmin, idcs.idcs_z[9]] = -1*Float64.(kron(I(N), I(idcs.n_sx)))
    b[idcs.idcs_bxmin] = Float64.(kron(ones(N), scppbm.xLowThd))


    A[idcs.idcs_bumax, idcs.idcs_z[2]] = Float64.(kron(I(N), scppbm.I_ul))
    A[idcs.idcs_bumax, idcs.idcs_z[10]] = Float64.(kron(I(N), I(idcs.n_su)))
    b[idcs.idcs_bumax] = Float64.(kron(ones(N), scppbm.uHighThd))
    A[idcs.idcs_bumin, idcs.idcs_z[2]] = Float64.(kron(I(N), scppbm.I_ul))
    A[idcs.idcs_bumin, idcs.idcs_z[11]] = -1*Float64.(kron(I(N), I(idcs.n_su)))
    b[idcs.idcs_bumin] = Float64.(kron(ones(N), scppbm.uLowThd))

    A[idcs.idcs_bpmax, idcs.idcs_z[3]] = Float64.(I(idcs.num_sp))
    A[idcs.idcs_bpmax, idcs.idcs_z[12]] = Float64.(I(idcs.num_sp))
    b[idcs.idcs_bpmax] = Float64.(scppbm.pHighThd)
    A[idcs.idcs_bpmin, idcs.idcs_z[3]] = Float64.(I(idcs.num_sp))
    A[idcs.idcs_bpmin, idcs.idcs_z[13]] = -1*Float64.(I(idcs.num_sp))
    b[idcs.idcs_bpmin] = Float64.(scppbm.pLowThd)

end 

function scp_init_l1cost!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing      # initial and parse L1-norm cost
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    N = scppbm.scpPrs.N
    idcs = subpbm.idcsLnrConPgm
    A, b = subpbm.A, subpbm.b

    # set the (idx_bxcn_bxcp, idx_x)
    A[idcs.idcs_bxcn, idcs.idc_z[1]] = Float64.(kron(I(N), scppbm.I_xc))
    A[idcs.idcs_bxcp, idcs.idc_z[1]] = Float64.(kron(I(N), scppbm.I_xc))

    # set the (idx_bxcn, idx_zsvc1_zsvc2)   theta + sxc1 - sxc2 = 0
    A[idcs.idcs_bxcn, idcs.idc_z[14]] = Float64.(kron(I(N), I(idcs.num_xc)))
    A[idcs.idcs_bxcn, idcs.idc_z[15]] = -1*Float64.(kron(I(N), I(idcs.num_xc)))
    # set the (idx_bxcn, idx_zsvc1_zsvc2)   theta - sxc1 + sxc3 = 0
    A[idcs.idcs_bxcp, idcs.idc_z[14]] = -1*Float64.(kron(I(N), I(idcs.num_xc)))
    A[idcs.idcs_bxcp, idcs.idc_z[16]] = Float64.(kron(I(N), I(idcs.num_xc)))

    # set the b[idx_bxcn_bxcp]
    b[idcs.idcs_bxcn] = zeros(idcs.num_xc)
    b[idcs.idcs_bxcp] = zeros(idcs.num_xc)

end

function  scp_init_tr!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing          # initial and parse trust region
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    idcs = subpbm.idcsLnrConPgm
    A = subpbm.A
    b = subpbm.b
    xref, uref, pref = scppbm.xref, scppbm.uref, scppbm.pref
    idcs_z = idcs.idcs_z
    
    # parameters auxiliary variables
    A[idcs.idx_bptrn, idcs.idx_zp] = Float64.(I(idcs.dims_p))
    A[idcs.idx_bptrn, idcs.idx_zptr[1]] = -1*Float64.(I(idcs.dims_p))           #p-etap+mup1 = pref
    A[idcs.idx_bptrn, idcs.idx_zptr[2]] = Float64.(I(idcs.dims_p))
    b[idcs.idx_bptrn] = pref

    A[idcs.idx_bptrp, idcs.idx_zp] = Float64.(I(idcs.dims_p))
    A[idcs.idx_bptrp, idcs.idx_zptr[1]] = Float64.(I(idcs.dims_p))           #p-etap+mup1 = pref
    A[idcs.idx_bptrp, idcs.idx_zptr[3]] = -1*Float64.(I(idcs.dims_p))
    b[idcs.idx_bptrp] = pref

    # state trust region
    A[idcs.idx_bxtr, idcs_z[1]] = Float64.(kron(I(N), I(nx)))
    A[idcs.idx_bxtr, idcs_z[18]] = Float64.(kron(I(N), [zeros(nx,1) -1*I(nx)]))
    b[idcs.idx_bxtr] = [x for xk in xref for x in xk]

    # control trust region
    A[idcs.idx_butr, idcs_z[2]] = Float64.(kron(I(N), I(nu)))
    A[idcs.idx_butr, idcs_z[19]] = Float64.(kron(I(N), [zeros(nu,1) -1*I(nu)]))
    b[idcs.idx_butr] = [u for uk in uref for u in uk]

end   

function scp_init_pgm!(subpbm::ScpSubPbm, scppbm::ScpPbm)::Nothing
    nx, nu, np = trjpbm.dynmdl.nx, trjpbm.dynmdl.nu, trjpbm.dynmdl.np
    idcs = subpbm.idcsLnrConPgm
    idcs_z = idcs.idcs_z
    idcs_b = idcs.idcs_b 
    A = subpbm.A
    Achk = subpbm.Achk
    Asp = subpbm.Asp
    num_ic, num_fc, = idcs.num_ic, idcs.num_fc
    pgm = subpbm.pgmLnrCon

    # check A
    Achk = copy(A)
    Achk[isnan.(A)] .= 0.0
    #= set sparse A's architecture
    
    [   blj_ic  | blkvc |   |
        blk_dyn |       | 0 | 0 ]
        blk_fc  |       |   |
     ----------------------------   
                A[]
     ]
    =#
    blk_Aic =  sparse( [[fill(NaN, num_ic, nx)  zeros(num_ic, nx*(N-1))]  zeros(num_ic, idcs.dims_u)  fill(NaN, num_ic, idcs.dims_p)] )
    blk_Adyn = sparse( [fill(NaN, idcs.num_dyn, idcs.dims_x+idcs.dims_u+idcs.dims_p)] )
    blk_Afc =  sparse( [[zeros(num_fc, nx*(N-1)) fill(NaN, num_fc, nx)]   zeros(num_fc, idcs.dims_u)  fill(NaN, num_fc, idcs.dims_p)] )
    blk_vc =   sparse( Achk[idcs_b[1][1]:idcs_b[3][end], idcs_z[4]] )
    blk_zero = spzeros(Float64, length(1:idx_bfc[end]), length(idcs_z[5][1]:idcs_z[end][end]) )

    blk_K = sparse( Achk[idx_bvcn[1]:end, :] )

    Asp = [ [ blk_Aic;
              blk_Adyn;
              blk_Afc;]   blk_vc  blk_zero ;

                        blk_K               ]
    
    subpgm.I, subpgm.J, _ = findnz(Asp)
    #for i  in 1:length(Asp.nzval)
    #    Asp[i] = A[I[i], J[i]]
    #end
   
    # set G
    dims_K = idcs.dims_K0 + idcs.dims_K2
    dims_v1 = idcs.dims_x+idcs.dims_u+idcs.dims_p+idcs.dims_vcdyn
    subpbm.G = [zeros(dims_K, dims_v1)  Float64.(I(dims_K))]
    subpgm.Gsp = sparse(G)

end

