"""
Overtake trajectory planning case-problem class
"""


# Main and Public functions
function TrjPbm_Ovrtk_Inin(

)
""" I. Overtake formula """
""" 1.1 trjPbm formula""" 
    """ A.1 Update Road&Lane constraints from HD Map, SLAM and Perception """
    #Local coordinate system repositioning
    """ B.1 Update Static obstacles from Perception"""
    """ B.2 Update Moving obstacles prediction and initial guess from perception and game layer """
    """ A.2 Update env and road parameters along guess from HD and VMC"""
    """ B.3 Update Dynamic parameters and constraints"""
    vehicledynamicmgr.Upd(mdl, constraints);
""" 1.2 OPC_SCP Pre-Hand-Parser""" 
    sc



end

function TrjPbm_Ovrtk_Main(

)

""" II. Long-peroidic/Events reformula in 1-5s or 10-30m 
              with Road&Lane, Env&Road update """
    # Event-triggered reformula: 1) driving state transition 2) trjPbm structure changes:
    # a.  
    # 1.Except Completely different driving scenarios, from highway to ramp merging  
    # trjPbmbe shall be incrementally constructed if there is any structure changes;   
    # 2.  shall be conA large set of Consecutive driving scenariosstructed into 
    # one flexible and redundant trjPbm, such as highway {acc, overtake, aeb};
    """ 1.1 trjPbm formula""" 
    # trjPbm and SCP structure are fixed. 
    # Only constraints info and some paras are updated from perception and VMC;

    # Construct the definition and data-structure
    trjPbmOvtk = TrjPbm_Ovrtk_Inin();
    # Select SCP algorithm and convex subproblem
    #if trjPbmOvtk.scpAgm == PTR
        scpParas = ScpParas(trjPbmOvtk);   
    #if trjPbmOvtk.cvxSubPbmType == lscop
        subPbm = ScpSubPbm(scpParas, trjPbmOvtk);

    #Other default parameters

    #Other data-structure
    soluscp = ScpSolu()
    histScp = ScpHist()

    """ 1.2 OPC_SCP Pre-Hand-Parser""" 
    # Initialize and pre-parse all problem and program


""" III. MPC: short-horizon periodic formula and solve in 100ms 
              with collision-free, game layer, auto dyn paras """

    """ 3.1 trjPbm formula of running perception""" 
        """ B.1 Update Static obstacles from Perception"""
        """ B.2 Update Moving obstacles prediction and initial guess from perception and game layer """
        """ B.3 Update Dynamic parameters and constraints"""
        TrjPbm_Upd_StcObs!(trjPbmOvtk);
        TrjPbm_Upd_MvgObs!(trjPbmOvtk);


    """ 3.2 OPC_SCP Hand-Parser of running perception"""
    



    """ 3.3 SCP Solver Loop: ScpSubProblem creating and solving in 10ms"""
        scpsubpbm_solve!(subPbm);
        scp_upd_trj!(subPbm, trjPbmOvtk);

        scp_recorder!(subPbm, trjPbmOvtk);

end






# Initialization and private functions

# I. Define the vehicle model and parameters
# Model A: simplest DubinCar







# end: Vehicle model and parameters 


# II. Collision-free


