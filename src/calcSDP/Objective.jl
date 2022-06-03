# Calculating the objective vector 
include("Cmatrixsigmatau.jl")
include("CalcProduct.jl")

# Sven implementation:

# cycles = collect(permutations(Array(Int8.(2:m))))
# for cycle in cycles
#     pushfirst!(cycle, 1)
# end

# distanceDict = Csigmatau(m)
# CmatrixDict = Dict() 
# # make distances from vertex to inverse of destination
# for (vertex, distance) in distanceDict 
#     vertexnew = deepcopy(vertex)
#     vertexnew[2:end] = vertexnew[end:-1:2];
#     CmatrixDict[vertexnew] = distance
# end

# NumbersToCycleCoefficients = Dict()
# for cycle in cycles 
#         # we compute C((1,...,m), cycle). 
#     currentOrbitNumber = labelCanonical(cycle)
#     NumbersToCycleCoefficients[currentOrbitNumber] = get!(NumbersToCycleCoefficients, currentOrbitNumber, 0) + (CmatrixDict[cycle])
# end

# display(NumbersToCycleCoefficients)

# New implementation:

# distance between inverse(cycle) and 1:n
# function distance(cycle)
#     c = [1, cycle[end:-1:2]]
# end

distanceDict = Csigmatau(m)
CmatrixDict = Dict() 
# make distances from vertex to inverse of destination
for (vertex, distance) in distanceDict 
    vertexnew = deepcopy(vertex)
    vertexnew[2:end] = vertexnew[end:-1:2];
    CmatrixDict[vertexnew] = distance
end

objective(cycle) = orbitSize(cycle) * CmatrixDict[cycle]

# Dict([cycle => objective(cycle) for cycle in keys(NumbersToCycleCoefficients)])

##
# (c,d) = rand(distanceDict)

# n = length(c)
# M = BitArray(undef, n,n)
# M .= 0
# for i = 1:n
#     M[i, c[i]] = 1
# end

# display(M')
# display(d)
# # display(c)

# display([(c[((i%m == 0) ? m : i%m)] - c[i-1] + m)%m - 1 for i = 2:m+1])
