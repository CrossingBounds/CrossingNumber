# Bounds for the crossing number of the complete bipartite graph

Computes the bounds $\alpha_m$ and $\beta_m$ for the crossing number of $K_{n,m}$ described in the paper

> *New lower bounds on crossing numbers of Km,n from
permutation modules and semidefinite programming*, Daniel Brosch and Sven Polak, ArXiv:?????, 2022.

To compute the bounds: 
- [SolveSDP_Alpha.jl](src/SolveSDP_Alpha.jl) contains a function `solveAlpha(m)` which computes $\alpha_m$.
- [SolveSDP_Beta_Iterative.jl](src/SolveSDP_Beta_Iterative.jl) contains a function `solveBetaIterative(m)` which computes $\beta_m$.
- [VisualizeSolutions.jl](src/VisualizeSolutions.jl) computes $\beta_m$ for $m\in\{5,7,9,11,13\}$ and plots the coefficients of the eigenvector of the rank 1 solutions.
