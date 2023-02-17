# Computes alpha_m with sdpa_dd.

include("calcSDP/BuildSDP.jl")

using SDPAFamily, LinearAlgebra, GenericLinearAlgebra
# using Hypatia, LinearAlgebra, GenericLinearAlgebra

function roundRationalPSD(A)
    if size(A,1) == 1
        return max(Rational(A), Rational{BigInt}(0//1))
    end
    eg = eigen(Hermitian(A))
    egVals = Rational.(eg.values)
    egVecs = Rational.(eg.vectors)
    posInds = egVals .> 0 
    return egVecs[:, posInds] * diagm(egVals[posInds]) * egVecs[:, posInds]'
end

"""
    solveAlpha(m)

Computes the bound alpha_m.
"""
function solveAlpha(m; variant = :sdpa_gmp)

    # setprecision(128)
    setprecision(512)
    
    @time sdpData = calcSDP(m);
    
    @time (t, yD, cD) = buildDualSDP_GMP(sdpData);
    @show length(cD)
    @show length(yD)

    
    
    constraints = [c[2] >= 0 for c in cD]
    P = maximize(t, constraints...; numeric_type=BigFloat)

    if variant == :sdpa_gmp
        params = SDPAFamily.Params{:sdpa_gmp,BigFloat}(
            maxIteration=200,
            epsilonStar=1e-40,
            lambdaStar=1e5,
            omegaStar=2,
            lowerBound=-1e5,
            upperBound=1e5,
            betaStar=0.1,
            betaBar=0.2,
            gammaStar=0.9,
            epsilonDash=1e-40,
            precision=512,
        )
        @time solve!(P, () -> SDPAFamily.Optimizer(presolve=false, params=params, verbose=SDPAFamily.VERBOSE; variant = :sdpa_gmp))

    elseif variant == :sdpa_dd 

        params = SDPAFamily.Params{:sdpa_dd,BigFloat}(
            maxIteration =    1000,
            epsilonStar  = 1.000e-17,
            lambdaStar   = 1.000e+03,
            omegaStar    = 2.000e+00,
            lowerBound   = -1.000e+05,
            upperBound   = 1.000e+05,
            betaStar     = 1.000e-01,
            betaBar      = 2.000e-01,
            gammaStar    = 9.000e-01,
            epsilonDash  = 1.000e-17
        )
    
        @time solve!(P, () -> SDPAFamily.Optimizer(presolve=false, params=params, verbose=SDPAFamily.VERBOSE; variant = :sdpa_dd))
        # @time solve!(P, () -> SDPAFamily.Optimizer())
    else
        @error "Use sdpa gmp or dd"
    end
    
    # @time solve!(P, () -> Hypatia.Optimizer{BigFloat}())


    @show P.optval
    @show P.status

 

    @show 8 * P.optval / (m * (m - 1))

    # rounding to a rational feasible solution
    # guaranteed PSD
    ySol = Dict(mu=>roundRationalPSD(Convex.evaluate(y)) for (mu, y) in yD)

    # potentially slightly too big
    tSol = Rational(Convex.evaluate(t))

    # checking constraints 
    @show tSol
    for v in sdpData.vars
        vVec = [i for i in v]
        A = -sum(dot(Symmetric(sdpData.coeff[mu][vVec]), ySol[mu]) for mu in keys(sdpData.coeff) if haskey(sdpData.coeff[mu],vVec); init = BigInt(0)//BigInt(1)) + sdpData.obj[v]
        B = sdpData.orbitSizes[v]
        # we need: A >= tSol*B
        # <=> tSol <= A/B
        @show tSol <= A//B
        tSol = min(tSol, A//B)   
    end
    @show tSol, Convex.evaluate(t)
    @show BigFloat(tSol) - Convex.evaluate(t)
    


    # dual 
    # @show sum(tr(Convex.evaluate(y)) for y in values(yD))

    # for (mu, blks) in sdpData.coeff

    #     @show mu
    #     Amu = sum(Symmetric(blks[Int8[j for j in c]])*P.constraints[i].dual for (i, c) in enumerate(sdpData.vars) if haskey(blks, Int8[j for j in c]))
    #     @show Amu 
    #     @show minimum(eigen(Hermitian(Amu)).values)
    # end

    
    # return P
    # return yD
    # return P.optval
    return tSol
end
