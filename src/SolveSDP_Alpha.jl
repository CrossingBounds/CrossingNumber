# Computes alpha_m with sdpa_gmp.

include("calcSDP/BuildSDP.jl")

using SDPAFamily, LinearAlgebra, GenericLinearAlgebra
# using Hypatia, LinearAlgebra, GenericLinearAlgebra

function roundRationalPSD(A)
    if size(A,1) == 1
        return max(Rational(A[1,1]), Rational{BigInt}(0//1))
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
    setprecision(128)
    
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

function solveAlphaDual(m; variant = :sdpa_gmp)
    setprecision(128)
    
    @time sdpData = calcSDP(m);
    
    @time (P, y, cD) = buildPrimalSDP_GMP(sdpData);
    # @show length(cD)
    # @show length(yD)

    
    
    # constraints = [c[2] >= 0 for c in cD]

    @info "Calling the solver..."

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
        @time solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=params, verbose=SDPAFamily.VERBOSE; variant = :sdpa_gmp))

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
    
        @time solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=params, verbose=SDPAFamily.VERBOSE; variant = :sdpa_dd))
        # @time solve!(P, () -> SDPAFamily.Optimizer())
    else
        @error "Use sdpa gmp or dd"
    end
    
    # @time solve!(P, () -> Hypatia.Optimizer{BigFloat}())


    @show P.optval
    @show P.status

 

    @show 8 * P.optval / (m * (m - 1))

    # # rounding to a rational feasible solution
    # # guaranteed PSD
    # ySol = Dict(mu=>roundRationalPSD(Convex.evaluate(y)) for (mu, y) in yD)

    # # potentially slightly too big
    # tSol = Rational(Convex.evaluate(t))

    # # checking constraints 
    # @show tSol
    # for v in sdpData.vars
    #     vVec = [i for i in v]
    #     A = -sum(dot(Symmetric(sdpData.coeff[mu][vVec]), ySol[mu]) for mu in keys(sdpData.coeff) if haskey(sdpData.coeff[mu],vVec); init = BigInt(0)//BigInt(1)) + sdpData.obj[v]
    #     B = sdpData.orbitSizes[v]
    #     # we need: A >= tSol*B
    #     # <=> tSol <= A/B
    #     @show tSol <= A//B
    #     tSol = min(tSol, A//B)   
    # end
    # @show tSol, Convex.evaluate(t)
    # @show BigFloat(tSol) - Convex.evaluate(t)
    


    # # dual 
    # # @show sum(tr(Convex.evaluate(y)) for y in values(yD))

    # # for (mu, blks) in sdpData.coeff

    # #     @show mu
    # #     Amu = sum(Symmetric(blks[Int8[j for j in c]])*P.constraints[i].dual for (i, c) in enumerate(sdpData.vars) if haskey(blks, Int8[j for j in c]))
    # #     @show Amu 
    # #     @show minimum(eigen(Hermitian(Amu)).values)
    # # end

    
    # # return P
    # # return yD
    # # return P.optval
    # return tSol
end

function writeSDPA(io::IO, m; precision = 256, sdpData = calcSDP(m))
    setprecision(precision)

    # sdpData = calcSDP(m)

    n = length(sdpData.vars)

    # last variable is linearly dependent
    vars = copy(sdpData.vars)[1:end-1]
    
    # last var = scaled 1 - other vars 
    finalOrbitSize = BigInt(sdpData.orbitSizes[sdpData.vars[end]])
    # finalVar = 1/finalOrbitSize - sum(y[i]*sdpData.orbitSizes[v] for (i,v) in enumerate(sdpData.vars[1:n-1]))/finalOrbitSize

    # x_i >= 0
    varBlocks = [:vars => -length(sdpData.vars)]

    # block (m) is redundant
    sdpBlocks = [mu => size(rand(b)[2],1) for (mu, b) in sdpData.coeff if mu != (AbstractAlgebra.Partition([length(sdpData.vars[1])]),1)]
    blockSizes = Dict(sdpBlocks)

    blocks = vcat(varBlocks, sdpBlocks)

    blockData = Dict(mu => Dict(c => [] for c in vcat(vars, [:const])) for (mu,) in blocks)
    for (i,c) in enumerate(vars)
        push!(blockData[:vars][c], (i, i, 1))
        push!(blockData[:vars][c], (n, n, -sdpData.orbitSizes[c]//finalOrbitSize))
    end
    push!(blockData[:vars][:const], (n,n, -1//finalOrbitSize))

    cEnd = [i for i in sdpData.vars[end]]
    for (mu, b) in sdpData.coeff
        mu == (AbstractAlgebra.Partition([length(sdpData.vars[1])]),1) && continue
        blk = zeros(Rational{Int}, blockSizes[mu], blockSizes[mu])
        for c in vars
            cv = [i for i in c]

            blk .= 0
            if haskey(b, cv)
                blk .= b[cv]
            end

            if haskey(b, cEnd)
                blk .-= 1//finalOrbitSize *sdpData.orbitSizes[c]*b[cEnd]
            end
            # @show blk

            for nonZ in findall(blk .!= 0)
                # @show nonZ 
                # @show blk 
                # @show blk[nonZ]
                push!(blockData[mu][c], (nonZ[1], nonZ[2], blk[nonZ]))#nonZ[1] == nonZ[2] ? blk[nonZ] : blk[nonZ]//2))
            end
        end
        if haskey(b, cEnd)
            blk = b[cEnd]
            for nonZ in findall(blk .!= 0)
                # @show nonZ 
                # @show blk 
                # @show blk[nonZ]
                push!(blockData[mu][:const], (nonZ[1], nonZ[2], -blk[nonZ]*1//finalOrbitSize))
            end
        end
    end

    objective = Rational{Int}.([sdpData.obj[v] for v in sdpData.vars[1:end-1]])
    for (i,c) in enumerate(vars)
        objective[i] -= 1//finalOrbitSize *sdpData.orbitSizes[c]*sdpData.obj[sdpData.vars[end]]
    end

    objShift = 1//finalOrbitSize*sdpData.obj[sdpData.vars[end]]

    @info "writing sdpa file"
    println(io, "\"alpha_$m")
    println(io, "\"objective shifted by $(-objShift)")
    println(io, "\"variable order:")
    for v in vars 
        println(io, "\"$v")
    end
    println(io, "\"last variable $(sdpData.vars[end]) eliminated")
    println(io, "\"block order:")
    for b in collect(keys(blockData))
        println(io, "\"$b")
    end
    println(io, length(vars))
    println(io, length(blockData))
    for b in blocks 
        print(io, b[2])
        print(io, " ")
    end
    println(io)
    for c in objective
        print(io, Int(c))
        print(io, " ")
    end
    println(io)
    for (mu, b) in blockData
        k = findfirst(x->x[1]==mu, blocks)
        for (v, entries) in b
            vi = findfirst(x->x==v, vars)
            if vi === nothing 
                vi = 0 
            end
            for (i, j, c) in entries 
                @assert i<=j
                # @show c 
                # @show typeof(c)
                @assert typeof(c) in [Rational{BigInt}, Int, Rational{Int}]
                println(io, "$vi $k $i $j $(BigFloat(Rational{BigInt}(c)))")
            end
        end
    end
    @info "finished writing file"
    # println(io, [b[2] for b in blocks])

    return vars, blocks, objective, blockData, objShift
end

##
using JLD2

fullPath = SDPAFamily.Optimizer().binary_path
fullPath = fullPath[1:end-3]
fullPath *= "dd"

m = 10
sdpData = calcSDP(m)

open("data$m.dat-s", "w") do io 
    global data = writeSDPA(io, m; sdpData = sdpData)
end

@save "sdpData$m.jld2" sdpData data
# @load "sdpData$m.jld2" sdpData data

@info "starting sdpa gmp"
# m = SDPAFamily.Optimizer()
# SDPAFamily.sdpa_gmp_binary_solve!(m, "data.dat-s", "out.dat")

@time run(pipeline(`$fullPath data$m.dat-s out$m.dat -pt 2`, stdout=stdout, stderr=stdout), wait=true)

output = readlines("out$m.dat")

objLine = output[findfirst(startswith("objValDual"), output)]
obj = parse(BigFloat,split(objLine, "= ")[2]) + data[end]

yMatLines = output[findfirst(startswith("yMat"), output)+2:end-4]
# yMatText = join(yMatLines)

curLine = 1

setprecision(256)

ySol = Dict()
for (mu, b) in data[2]
    if b < 0
        yMatText = yMatLines[curLine]
        curLine += 1
    else
        yMatText = join(yMatLines[curLine:curLine+b-1])
        curLine += b 
    end
    yMatText = replace(yMatText, "{"=>"", "}"=>"", " "=>"")
    if yMatText[end] == ','
        yMatText = yMatText[1:end-1]
    end

    # @show yMatText
    ySplit = split(yMatText, ",")
    yDat = parse.(BigFloat, ySplit)
    if b > 0
        yDat = reshape(yDat, b,b)
    end
    ySol[mu] = yDat
    # @show mu 
    # display(yDat)
end

# rounding to a rational feasible solution
# guaranteed PSD
ySol = Dict(mu=>roundRationalPSD(y) for (mu, y) in ySol if mu != :vars)
ySol[(AbstractAlgebra.Partition([m]),1)] = [0//1;;]

# potentially slightly too big
tSol = Rational(obj)

setprecision(512)

# checking constraints 
@show tSol
minGap = 1000
for v in sdpData.vars
    vVec = [i for i in v]
    # -sum(dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), Y[mu]) for mu in keys(sdpData.coeff) if mu in blocks && haskey(sdpData.coeff[mu],[i for i in v]);init=0) - t*sdpData.orbitSizes[v] + sdpData.obj[v] >= 0
    A = -sum(dot(Symmetric(sdpData.coeff[mu][vVec]), ySol[mu]) for mu in keys(sdpData.coeff) if haskey(sdpData.coeff[mu],vVec); init = BigInt(0)//BigInt(1)) + sdpData.obj[v]
    B = sdpData.orbitSizes[v]
    # we need: A >= tSol*B
    # <=> tSol <= A/B
    # @show tSol <= A//B
    minGap = min(minGap, A//B - tSol)
    tSol = min(tSol, A//B)   
end
@show tSol, obj
@show BigFloat(tSol) - obj
@show BigFloat(minGap)