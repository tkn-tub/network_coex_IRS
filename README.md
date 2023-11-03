# network_coex_IRS

This is built upon [Matlab - Reconfigurable Intelligent Surfaces (RIS)](https://de.mathworks.com/help/phased/ug/introduction-to-reconfigurable-intelligent-surfaces.html) example, following the 
 proposed model presented in [1].

While in [1] the environment is represented in the 2-D space, the example allows defining a 3-D environment. Multiple parameters can be defined, such as the number of operators, subsurfaces, elements of IRS subsurfaces, random trials, number of drops etc.

This example has been extended to:
- 'Random Assignments' - which computes the sum rate of the whole network of operators randomly assigned to a subsurface.
 

## System specification:
'helperRISSurface' is used in the example, thus the 'IRS_sub_rand.m' file should be located on the same path with 'MATLAB/Examples/R2023b/phased/IntroductionToRISExample'.

## References
[1] J. Angjo, A. Zubow and F. Dressler, "Coexistence Challenges in IRS-assisted Multi-Operator Networks," Proceedings of IEEE Global Telecommunications Conference (GLOBECOM 2023), 4th Workshop on Emerging Topics in 6G Communications (6GComm), Kuala Lumpur, Malaysia, December 2023. (to appear) \
