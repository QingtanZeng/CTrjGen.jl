Julia implemetation of PTR for vehicle trajectory generation, solely for academic purposes,
which take 8-Dof Auto elevated-road cases and 3-Dof rocket landing cases with hand-parsed process.

Abbreviation
PTR: 

Implementation Highlights
1. Hand-Parsed process into
  a) linear conic programming using second-order interior-points method (homogenous self-dual embedding IPM, ECOS)
  b) generic conic programming using first-order primal-dual method, including SCS(ADMM+IPM) and PIPG

2. Model common objective and constraints from planning and control of vehicle
  1) Auto high-speed cases:
    a.  
  3) Rocket landing cases: 
