

abstract type AbstTrjOpc end

mutable struct AutoTrjOpc <: AbstTrjOpc

    # 2.0 Constraints
    # 2.1 Dynami system
    dynmdl::DynMdl
    # 2.2, 2.3 states' Ampl-Box and control Ampl&Slop constraints without extra need
    dyncstr::DynCstr
    #2.3 boundaries
    initbrdy::DynBdry
    termbrdy::DynBdry

    #2.4 extra states constraits
    # 2.4.1 collision-free states constraints


end