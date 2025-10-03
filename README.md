the julia implemetation of PTR Algorithms for vehicle trajectory generation, which take 8Dof Auto high-speed cases and 3-Dof rocket landing cases.

1. Hand-Parsed process with
  1) linear conic programming using second-order interior-points method (homogenous self-dual embedding IPM, ECOS)
  2) generic conic programming using first-order primal-dual method, including SCS(ADMM+IPM) and PIPG

2. Model common objective and constraints from planning and control of vehicle
  1) Auto high-speed cases:
    a.  
  3) Rocket landing cases: 
