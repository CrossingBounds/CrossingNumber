## OLD VARIANT: READS SDPA FILES


using JuMP
using MosekTools
using MathOptInterface

# Solve primal and dual SDPs using JuMP

# m = read_from_file("/mnt/c/users/dbrosch/Downloads/sdp_crossing_m10_S10.dat-s")
set_optimizer(m, Mosek.Optimizer)
optimize!(m)

# list_of_constraint_types(m)
tmp = rand(all_constraints(m, Vector{AffExpr}, MathOptInterface.PositiveSemidefiniteConeTriangle))

##

# A = load("CrossingM7SymBasis.jld2")["basis"]
# B = load("CrossingM7SymBasisOld.jld2")["basis"]

## Solving problems calculated by Sven

using JuMP
using SparseArrays

function loadSparseSDPA(filename)
    numVars = 0
    nBlks = 0
    blkSizes = Int64[]
    obj = Int64[]
    blks = Dict()

    i = 1

    for ln in eachline(filename)
        # @show ln
        if length(ln) == 0
            continue
        end
        if i == 1
            numVars = parse(Int, ln)
        elseif i == 2
            nBlks = parse(Int, ln)
        elseif i == 3
            blkSizes = parse.(Int, split(ln, " ")[1:end-1])
        elseif i == 4
            obj = parse.(Int, split(ln, " ")[1:end-1])
        else
            matno, blkno, p, q, c = parse.(Int32, split(ln, " "))

            if !haskey(blks, matno)
                blks[matno] = Dict()
            end
            if !haskey(blks[matno], blkno)
                blks[matno][blkno] = spzeros(Int32, blkSizes[blkno], blkSizes[blkno])#[spzeros(Int32, blkSizes[i], blkSizes[i]) for i = 1:nBlks]
            end
            blks[matno][blkno][p,q] = c
        end

        i += 1
    end

    return (numVars = numVars, nBlks = nBlks, blkSizes = blkSizes, obj = obj, blks = blks)
end

# loadSparseSDPA("sdp_crossing_m4_S4xS2.dat-s")

function makeModelFromSDPA(filename)#, nonnegVars = true)
    data = loadSparseSDPA(filename)

    m = Model()

    # y = nonnegVars ? @variable(m, y[1:data.numVars] >= 0) : @variable(m, y[1:data.numVars])
    
    y = @variable(m, y[1:data.numVars] >= 0)

    @objective(m, Min, dot(y, data.obj))

    startBlk = data.numVars + 1 #nonnegVars ? data.numVars + 1 : 1

    for b = startBlk:data.nBlks
        @show (b, data.nBlks, data.blkSizes[b])

        #TODO: optimize to only calc triangle

        blk = Symmetric([ k<=l ? [haskey(data.blks[i],b) ? data.blks[i][b][k,l] : 0 for i = 1:data.numVars]' * y : 0*y[1] for k = 1:data.blkSizes[b], l = 1:data.blkSizes[b] ])

        # blks[mu] = Symmetric([ [haskey(m.sdpData[G],mu) ? m.sdpData[G][mu][i,j] : 0//1 for G in keys(m.sdpData)]' * y for i = 1:length(m.sdpBasis[mu]), j = 1:length(m.sdpBasis[mu])])
        
        # blk = sum(y[i] * Symmetric(data.blks[i][b]) for i = 1:data.numVars if haskey(data.blks[i],b))

        # blks[mu] = sum(y[G] * Symmetric(m.sdpData[G][mu]) for G in keys(m.sdpData) if haskey(m.sdpData[G],mu))



        if  haskey(data.blks[0],b)
            blk -=  Symmetric(data.blks[0][b])
        end

        if data.blkSizes[b] > 1
            @constraint(m, blk in PSDCone())
        else
            @constraint(m, blk .>= 0)
        end
    end

    return (m,y)
end

function makeModelFromSDPAPrimal(filename)
    data = loadSparseSDPA(filename)

    m = Model()

    Y = []

    for (i,b) in enumerate(data.blkSizes)
        @show (i, length(data.blkSizes))
        if b == 1
            push!(Y, @variable(m, base_name = "Y$i"))
            @constraint(m, Y[i] >= 0)
        else
            push!(Y, @variable(m, [1:b, 1:b] in PSDCone(), base_name = "Y$i"))
            # @constraint(m, Y[i] in PSDCone())
            # @constraint(m, Y[i] .== 0)
        end
    end

    @objective(m, Max, sum(dot(Symmetric(data.blks[0][i]), Y[i]) for i = 1:data.nBlks if haskey(data.blks[0],i)))

    for i = 1:data.numVars
        @show (i, data.numVars)
        @constraint(m, sum(dot(Symmetric(data.blks[i][j]), Y[j]) for j = 1:data.nBlks if haskey(data.blks[i],j)) == data.obj[i])
    end

    return (m,Y)
end

# m = makeModelFromSDPA("sdp_crossing_m4_S4xS2.dat-s")

##

using LinearAlgebra

k = 9

# blks = loadSparseSDPA("sdp_crossing_m$(k)_S$(k)xS2.dat-s")

# 77579220 (sparse Int32)
# 113870336 (full Int32)
# Base.summarysize(blks.blks)

#
@time (m,y) = makeModelFromSDPAPrimal("data/sdp_crossingWindows_m9_S9xS2.dat-s")
# @time (m,y) = makeModelFromSDPAPrimal("sdp_crossing_m$(k)_S$(k)xS2.dat-s")
# @time (m,y) = makeModelFromSDPA("sdp_crossing_m$(k)_S$(k).dat-s")

# using JLD2
# using FileIO

# save("ModelM$(k)Dual.jld2","my",(m,y))

##

# (m,y) = load("ModelM10Dual.jld2")["my"]

using MosekTools
set_optimizer(m, Mosek.Optimizer)

@time optimize!(m)
@show objective_value(m)
@show primal_status(m)
@show dual_status(m)
@show termination_status(m)

display(minimum([minimum(value.(x)) for x in y]))
display(maximum([maximum(value.(x)) for x in y]))

display(sum(tr(value.(x)) for x in y))


# for x in y
#     if maximum(abs.(value.(x))) > 0.00001
#         display(x)
#         # display(value.(x))
#     end
# end

display(sum(sum(abs.(value.(x)) .> 0.00001) for x in y))

# display(sum(sum(abs.(value.(x)) .> 0.001) for x in y[1:239]))
display(sum(typeof(x) == VariableRef ? 1 : length(x) for x in y))

## SDPNAL test: Needs upper and lower variable bounds, is slow...
using SDPNAL
using MATLAB
mat"path(pathdef)"

@time (m,y) = makeModelFromSDPAPrimal("sdp_crossing_m$(k)_S$(k)xS2.dat-s", true)
set_optimizer(m, SDPNAL.Optimizer)
# set_optimizer_attribute(m, "maxiter", 1000)
# set_optimizer_attribute(m, "printlevel", 2)

# for x in y
#     @constraint(m, x .<= 140)
#     @constraint(m, x .>= -1)
# end

@constraint(m, sum(typeof(x) == VariableRef ? x : tr(x) for x in y) <= 45000)

println("Starting SDPNAL")

@time optimize!(m)

@show objective_value(m)
@show primal_status(m)
@show dual_status(m)
@show termination_status(m)



##
using JuMP, SDPNAL
m = Model(SDPNAL.Optimizer)
@variable(m, Z)
@constraint(m, Z >= 0)
@objective(m, Min, Z)
optimize!(m)