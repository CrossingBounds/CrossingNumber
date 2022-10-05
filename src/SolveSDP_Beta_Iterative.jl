include("calcSDP/BuildSDP.jl")


using SDPAFamily, LinearAlgebra#,JLD2

"""
    solveBetaIterative(m)

Computes the bound beta_m and solves it iteratively.
"""
function solveBetaIterative(m)

    
    mu = (Partition([m - 2, 1, 1]), -1)
    dualBlocksRestricted = [mu]


    # sdpData = load("data/sdpDataBeta$m.jld2")["sdpData"]

    @time sdpData = calcSDP(m, dualBlocksRestricted);
    # save("data/sdpDataBeta$m.jld2", "sdpData", sdpData)

    @time (t, yD, cD) = buildDualSDP_GMP(sdpData, dualBlocksRestricted);

    addedConstraints = Set()
    addedConstraintsDict = Dict()
    addedConstraintsInd = Set()
    bigM = 1e+5
    constraints = Any[t<=bigM]


    for m in keys(yD)
        push!(constraints, sum(yD[m]) >= -100)
    end

    it = 1
    maxIt = 100

    setprecision(256)

    @time while true
        if it == maxIt
            @warn "Reached iteration limit, aborting..."
            break
        end
        @info("Iteration $(it)")
        it += 1

        P = maximize(t, constraints...; numeric_type=BigFloat)

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

        @time solve!(P, () -> SDPAFamily.Optimizer(presolve=false, params=params))

        @show P.optval
        @show P.status

        global optVal = P.optval

        maxViolation = 1
        ind = 0
        lastY = Convex.evaluate(yD[mu])
        lastT = Convex.evaluate(t)

        
        @showprogress "checking constraints" for (i,v) in enumerate(sdpData.vars)
            if !(i in addedConstraintsInd)
                constr = 0
                if haskey(sdpData.coeff[mu],[i for i in v])
                    constr = -dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), lastY)
                end
                constr += - lastT*sdpData.orbitSizes[v] + sdpData.obj[v]

                if constr < maxViolation
                    maxViolation = constr 
                    ind = i
                end
            end
        end


        @show maxViolation
        if maxViolation >= 0.000001
            @info("SDP Solved!")
            break
        else
            @info("Adding constraint $ind")
            push!(addedConstraints, cD[ind])
            addedConstraintsDict[ind] = cD[ind][2] >= 0
            push!(constraints, cD[ind][2] >= 0)
            push!(addedConstraintsInd, ind)
        end
    end

    @show 8 * optVal / (m * (m - 1))

    return optVal, Symmetric(Convex.evaluate(yD[mu])), addedConstraintsDict
end
