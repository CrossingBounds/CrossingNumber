# Checking the decomposition conjecture for S_m symmetry alone

include("StuffFromFlagSOS.jl")
using LinearAlgebra, AbstractAlgebra

m = 4

lambda = Partition(repeat([1],m))
rowAut = (size = m, gen = [vcat(2:m,1)])

blkSizes = []
blks = []
bs = biggerShapes(lambda)
semistandards = []

sdpBasisFull = Dict()

for (i,mu) in enumerate(bs)
    K = Kostka(mu, lambda)

    reynolds = zeros(Int64, K[1], K[1])
    if rowAut.size > 1
        generator = []
        for p in rowAut.gen
            push!(generator, (p, semiTabPermMatV3(K, p)))
        end
        fullGroup = generateGroup(generator, rowAut.size, true)
        for p in fullGroup[1]
            reynolds += fullGroup[2][p]
        end
    else
        for i = 1:K[1]
            reynolds[i, i] = 1
        end
    end
    r = rank(reynolds)
    push!(blkSizes, r)
    push!(blks, reynolds)
    push!(semistandards, K[2])
end

# returns true if it should be in basis according to conjecture
function conjecture(t::AbstractAlgebra.Generic.YoungTableau{Int64})
    res = 0
    for i = 1:m-1
        if findfirst(x->x==i, t)[1] < findfirst(x->x==i+1, t)[1] 
            res += i
        end
    end

    return res % m == 0
end

function swapRole(t::AbstractAlgebra.Generic.YoungTableau{Int64})
    
    f = t.fill
    newfill = [findfirst(x->x==i, f) for i = 1:length(f)]
    @show f
    @show newfill
    return YoungTableau(t.part, newfill)
end

for i = 1:length(bs)
    A = blkSizes[i]
    chosen = conjecture.(semistandards[i])
    B = sum(chosen)
    C = rank(blks[i][:, chosen])

    
    display((bs[i],A,B,C))
    display.(semistandards[i][chosen])
    if A!=B
        @error("Number of tableaux wrong!")
    elseif A!=C 
        @error("Number correct, but linearly dependent!")
        display(blks[i])
        display(chosen)

    end
end

##


include("StuffFromFlagSOS.jl")

m = 11

lambda = Partition(repeat([1],m))
rowAut = (size = m, gen = [vcat(2:m,1)])

blkSizes = []
blks = []
bs = biggerShapes(lambda)
semistandards = []

sdpBasisFull = Dict()

for (i,mu) in enumerate(bs)
    K = Kostka(mu, lambda)
    push!(semistandards, K[2])
end

blkSizes = [sum(conjecture.(semistandards[i])) for i = 1:length(bs)]

@show sum(blkSizes)
@show sum(b^2 for b in blkSizes)