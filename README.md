# network_coex_IRS

This is built upon [Matlab - Reconfigurable Intelligent Surfaces (RIS)](https://de.mathworks.com/help/phased/ug/introduction-to-reconfigurable-intelligent-surfaces.html) example, following the 
 proposed model presented in [1].

While in [1] the environment is represented in the 2-D space, the example allows defining a 3-D environment. Multiple parameters can be defined, such as the number of operators, subsurfaces, elements of IRS subsurfaces, random trials, number of drops etc.

This example has been extended to:
- 'Random Assignments' - which computes the sum rate of the whole network of operators randomly assigned to a subsurface.
 

## System specification:
'helperRISSurface' is used in the example, thus the 'IRS_sub_rand.m' file should be located on the same path with 'MATLAB/Examples/R2023b/phased/IntroductionToRISExample'.

## References
[[1]](https://ieeexplore.ieee.org/document/10464696) J. Angjo, A. Zubow and F. Dressler, "Side Effects of IRS: On the Need for Coordination in 6G Multi-Operator IRS-Assisted Networks," 2023 IEEE Globecom Workshops (GC Wkshps), Kuala Lumpur, Malaysia, 2023, pp. 1380-1385, doi: 10.1109/GCWkshps58843.2023.10464696. \
