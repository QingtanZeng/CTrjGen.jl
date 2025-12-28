> :information_source: The repository CTrjGen is being developed in parallel with the author's R&D progress,
> to implement more new design features from academic papers from 2023-2026 and practical functions
> including AutoDiff, SCS solver and so on.

<p align="center">
<img alt="CTrjGen FlowChart"
    title="CTrjGen FlowChart"
    src="media/CTrjGen01.png"
    width="800px" />
</p>

<p align="center">
    <a href="http://www.gnu.org/licenses/gpl-3.0.txt"><img src="https://img.shields.io/badge/license-GPL_3-green.svg" alt="License GPL 3" /></a>
</p>

The <b>Computational Trajectory Generation</b> (CTrjGen) is a Julia implementation of the SCP PTR algorithm [1] 
for trajectory generation with a <b>hand-parsed process</b> and ECOS (a linear SOCP solver) [4], solely for academic purposes. 

Sequential Convex Programming (SCP) is a type of multiple-shooting direct method for numerical optimal control problems, 
modeled from autonomous systems such as Autonomous Driving, Robotic Loco-manipulation, Rocket Landing and so on. The Penalized 
Trust Region (PTR) is an SCP algorithm designed by Dr. T. P. Reynolds, Prof. B. Açıkmeşe et al. from ACL, University of Washington [1].

## Overview
CTrjGen is mainly composed of three parts.
1. OCP Class: dynamics, constraints, cost and ocp parameters.
2. SCP Class and parser function.
3. SubProblem class and lsocp structure used in ECOS.

## Implementation Highlights
1. <b>Hand-Parsed process</b>

- Rather than using automatic parsers such as JuMP or CVXPY [2], the standard form of linear SOCP is hand-parsed from the original OCP.
Because developers must be <b>familiar with each calculation step, eliminate computational bottlenecks and test software performance</b>,
- especially for real-time embedded systems, which usually
have limited computing resources, require functional safety review and failure may cause losses in the real world.

</p>
<p align="center">
<img alt="Sparse structure of A,b,c Array from linear SOCP"
    title="Sparse structure of A,b,c Array from linear SOCP"
    src="src/trajgen/examples/trjdb_pgm_init.png"
    width="400px" />
</p>

2. <b>Inverse-free FOH discretization</b> using RK4 of nonlinear system [3].

## Reference
[1] Reynolds, T. P. (2020). Computational guidance and control for aerospace systems. University of Washington.

[2] https://github.com/UW-ACL/SCPToolbox.jl

[3] Kamath, A. G., Doll, J. A., Elango, P., Kim, T., Mceowen, S., Yu, Y., ... & Açıkmeşe, B. (2025). Onboard Dual Quaternion Guidance for Rocket Landing. arXiv preprint arXiv:2508.10439.

[4] Domahidi, A., Chu, E., & Boyd, S. (2013, July). ECOS: An SOCP solver for embedded systems. In 2013 European control conference (ECC) (pp. 3071-3076). IEEE.

## License
Copyright (C) 2025 Qingtan Zeng

This project is licensed under the GPLv3 License - see the [LICENSE](LICENSE) file for details.
