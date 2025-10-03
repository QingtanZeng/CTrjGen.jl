Julia implemetation of PTR for vehicle trajectory generation, solely for academic purposes, which take 8-Dof Auto elevated-road cases and 3-Dof rocket landing cases with hand-parsed process.

Abbreviation
PTR:
SCP:
IPM:
PIPG:
ADMM:

Implementation Highlights
1. Hand-Parsed process into
  a) linear conic programming using second-order interior-points method (homogenous self-dual embedding IPM, ECOS)
  b) generic conic programming using first-order primal-dual method, including SCS(ADMM+IPM) and PIPG

2. Modelling common objective and constraints from planning and control of vehicle
  a) 8-Dof Auto elevated-road cases:
  
  b) 3-Dof rocket landing cases: 

References
[1] 
[2] 
[3] 
[4]