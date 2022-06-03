using JuMP
using FileIO

include("BlockDiagonalize.jl")
include("Cmatrixsigmatau.jl")
# include("Objective.jl")

function calcSDP(m, justBlocks = [])
    @info("Block diagonalize")
    @time basis = decomposeModule(m, justBlocks)
    @time coeff = blockDiagonalize(basis)
    @time vars = union!([Tuple.(keys(b)) for b in values(coeff)]...)

    @info("CSigmaTau")
    @time distanceDict = CsigmatauOrbits(m)
    CmatrixDict = Dict() 
    # make distances from vertex to inverse of destination
    for (vertex, distance) in distanceDict 
        vertexnew = deepcopy(vertex)
        vertexnew[2:end] = vertexnew[end:-1:2];
        CmatrixDict[labelCanonical(Int.(vertexnew))] = distance
    end
    @show length(vars)
    
    # for v in Tuple.(keys(distanceDict))
    #     if !(v in vars)
    #         @show v
    #         @show CmatrixDict[[i for i in v]]
    #     end
    # end
    
    vars = union!(vars, Tuple.(keys(distanceDict)))
    @show length(vars)

    orbitSizes = Dict(v=>orbitSize([i for i in v]) for v in vars)

    obj = Dict(v=>orbitSizes[v] * CmatrixDict[[i for i in v]] for v in vars)
    
    

    return (coeff= coeff, vars=Tuple.(vars),obj=obj, orbitSizes = orbitSizes)
end


# sdpData = calcSDP(4)

function buildPrimalSDP(sdpData)
    m = Model()

    y = @variable(m, y[sdpData.vars] >= 0)

    constraints = []

    push!(constraints, @constraint(m, sum(y[v]*sdpData.orbitSizes[v] for v in sdpData.vars) == 1))

    for (mu, b) in sdpData.coeff
        n = size(rand(b)[2],1)

        blk = Symmetric([ k<=l ? [haskey(b,[i for i in v]) ? b[[i for i in v]][k,l] : 0 for v in sdpData.vars]' * y : GenericAffExpr{Float64, VariableRef}() for k = 1:n, l = 1:n ])

        # psdBlocks = sum(blkD.blks[i] .* x[i] for i = 1:P.n)
        # N = size(blkD.blks[1],1)
        # psdBlocks = Symmetric([ k<=l ? [blkD.blks[i][k,l] for i in 1:P.n]' * x : GenericAffExpr{Float64, VariableRef}() for k = 1:N, l = 1:N ])

        if n == 1
            push!(constraints, (mu,@constraint(m, blk .>= 0)))
        else
            push!(constraints, (mu,@constraint(m, blk in PSDCone())))
        end
    end

    @objective(m, Min, sum(y[v] * sdpData.obj[v] for v in sdpData.vars))

    return (m, y, constraints)

end


"""
    buildDualSDP(sdpData, blocks = collect(keys(sdpData.coeff)), addConstraints = true)

Returns a jump model given by the sdpData. If blocks is a subset of partitions indexing the blocks, it generates the SDP with just the blocks corresponding to the given partitions. If addConstraints is true, it adds all (linear inequality) constraints to the SDP. If false, it still creates the constraints, and returns them, but does not add them to the SDP.
"""
function buildDualSDP(sdpData, blocks = collect(keys(sdpData.coeff)), addConstraints = true)
    m = Model()

    Y = Dict()
    
    constraints = []

    t = @variable(m, t)

    for (mu, b) in sdpData.coeff
        if mu in blocks
            n = size(rand(b)[2],1)
            if n == 1
                Y[mu] = @variable(m, base_name = "Y$mu")
                @constraint(m, Y[mu] >= 0)
            else
                Y[mu] = @variable(m, [1:n, 1:n] in PSDCone(), base_name = "Y$mu")
            end
        end
    end
    
    for v in sdpData.vars
        constrTerm = -sum(dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), Y[mu]) for mu in keys(sdpData.coeff) if mu in blocks && haskey(sdpData.coeff[mu],[i for i in v])) - t*sdpData.orbitSizes[v] + sdpData.obj[v]

        if addConstraints
            push!(constraints, (v,constrTerm, sdpData.obj[v], @constraint(m, constrTerm >= 0)))
        else
            push!(constraints, (v,constrTerm, sdpData.obj[v]))
        end
    end

    @objective(m, Max, t)

    return (m, t, Y, constraints)

end

##
using Convex
"""
    buildDualSDP_GMP(sdpData, blocks = collect(keys(sdpData.coeff)))

Returns a Convex model given by the sdpData. If blocks is a subset of partitions indexing the blocks, it generates the SDP with just the blocks corresponding to the given partitions..
"""
function buildDualSDP_GMP(sdpData, blocks = collect(keys(sdpData.coeff)))
    # m = Model()

    Y = Dict()
    
    constraints = []

    t = Variable(1)

    for (mu, b) in sdpData.coeff
        if mu in blocks
            n = size(rand(b)[2],1)
            if n == 1
                Y[mu] = Variable(1, Positive())
            else
                Y[mu] = Semidefinite(n)
            end
        end
    end
    
    for v in sdpData.vars
        constrTerm = -sum(Iterators.flatten((0,dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), Y[mu]) for mu in keys(sdpData.coeff) if mu in blocks && haskey(sdpData.coeff[mu],[i for i in v])))) - t*sdpData.orbitSizes[v] + sdpData.obj[v]

        push!(constraints, (v,constrTerm, sdpData.obj[v]))
    end

    # obj = 1*t#@objective(m, Max, t)

    return (t, Y, constraints)

end

# function buildRestrictedDualSDP(sdpData, blocks)
#     m = Model()

#     Y = Dict()
    
#     constraints = []

#     t = @variable(m, t)

#     for (mu, b) in sdpData.coeff
#         if mu in blocks
#             n = size(rand(b)[2],1)
#             if n == 1
#                 Y[mu] = @variable(m, base_name = "Y$mu")
#                 @constraint(m, Y[mu] >= 0)
#             else
#                 Y[mu] = @variable(m, [1:n, 1:n] in PSDCone(), base_name = "Y$mu")
#             end
#         end
#     end
    
#     for v in sdpData.vars
#         push!(constraints, (v,@constraint(m, sum(dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), Y[mu]) for mu in keys(sdpData.coeff) if mu in blocks && haskey(sdpData.coeff[mu],[i for i in v])) + t*sdpData.orbitSizes[v] <= sdpData.obj[v])))
#     end

#     @objective(m, Max, t)

#     return (m, Y, constraints)

# end