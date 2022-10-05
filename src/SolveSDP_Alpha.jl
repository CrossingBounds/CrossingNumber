# Computes alpha_m with sdpa_dd.

include("calcSDP/BuildSDP.jl")

using SDPAFamily, LinearAlgebra

"""
    solveAlpha(m)

Computes the bound alpha_m.
"""
function solveAlpha(m)

    setprecision(128)
    
    @time sdpData = calcSDP(m);
    
    @time (t, yD, cD) = buildDualSDP_GMP(sdpData);

    
    
    constraints = [c[2] >= 0 for c in cD]
    P = maximize(t, constraints...; numeric_type=BigFloat)

    params = SDPAFamily.Params{:sdpa_dd,BigFloat}(
        maxIteration=200,
        epsilonStar=1e-14,
        lambdaStar=1e5,
        omegaStar=2,
        lowerBound=-1e5,
        upperBound=1e5,
        betaStar=0.1,
        betaBar=0.2,
        gammaStar=0.9,
        epsilonDash=1e-14,
    )

    @time solve!(P, () -> SDPAFamily.Optimizer(presolve=false, params=params, verbose=SDPAFamily.VERBOSE))

    @show P.optval
    @show P.status

 

    @show 8 * optVal / (m * (m - 1))

    return optVal
end
