Julia implemetation of SCP PTR for trajectory generation, solely for academic purposes with hand-parsed process.

Abbreviation
SCP:
IPM:

Overview



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