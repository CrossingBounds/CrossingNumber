# Bounds for the crossing number of the complete bipartite graph

Computes the bounds $\alpha_m$ and $\beta_m$ for the crossing number of $K_{n,m}$ described in the paper

> *New lower bounds on crossing numbers of Km,n from
permutation modules and semidefinite programming*, Daniel Brosch and Sven Polak, [ArXiv:2206.02755](https://arxiv.org/abs/2206.02755), 2022.

To compute the bounds: 
- [SolveSDP_Alpha.jl](src/SolveSDP_Alpha.jl) contains a function `solveAlpha(m)` which computes $\alpha_m$ with [Convex.jl](https://github.com/jump-dev/Convex.jl) and rounds the solution to a feasible rational solution. The same file also provides a more efficient function `solveAlphaExternal(m)` which manually generates the SDPA file, avoiding Convex.jl and JuMP.jl, and verifies the solution.
- [SolveSDP_Beta_Iterative.jl](src/SolveSDP_Beta_Iterative.jl) contains a function `solveBetaIterative(m)` which computes $\beta_m$ and rounds the solution to a feasible rational solution.
- [VisualizeSolutions.jl](src/VisualizeSolutions.jl) computes $\beta_m$ for $m\in\{5,7,9,11,13\}$ and plots the coefficients of the eigenvector of the rank 1 solutions.

## Quick test of solveAlpha

Assume Julia installed; run `julia` in terminal,
in the main package directory (with `readme.md`).
- hit `]` to get `pkg>` prompt
- type `activate .` (and hit return)
- type `initialize`
- hit backspace to get back to `julia>` prompt
- type `include("src/SolveSDP_Alpha.jl")`
- type `solveAlpha(5)`; after a little while, it will give a rational close to `1.947213595...` as the answer.
- etc... (e.g., type `solveAlpha(7)` to get `4.359315494807...`)
