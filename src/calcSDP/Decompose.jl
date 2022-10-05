
# Decomposes the module into irreducible submodules
include("StuffFromFlagSOS.jl")



# New code for this package
using LinearAlgebra

# Full symmetry reduction
function decomposeModule(m, justBlocks = [])

    # if m==9
    #     return load("data/decompose9.jld2")["basis"]
    # end

    lambda = Partition(repeat([1],m))
    rowAut = (size = m, gen = [vcat(2:m,1)])

    blkSizes = []
    blks = []
    #allReynolds = Dict()
    # test = Dict([])
    bs = biggerShapes(lambda)

    if length(justBlocks) > 0
        bs = [b for b in bs if b in [a[1] for a in justBlocks]]
    end

    ## Full sym:
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
    end

    for (i,mu) in enumerate(bs)
        K = Kostka(mu, lambda)

        invMat = semiTabPermMatV3(K, collect(m:-1:1))
        R = blks[i]


        sym = R + invMat*R
        aSym = R - invMat*R

        #allReynolds[(mu, 1)] = sym
       # allReynolds[(mu, -1)] = aSym

        # test[mu] = (K, sym, aSym)

        # if rank(sym) + rank(aSym) > 0
        #     display((i,mu, rank(sym), rank(aSym)))
        # end

        inds = findColSpanningSetV3(sym)
        if length(inds) > 0

            sdpBasisFull[(mu,1)] = []

            for i in inds
                push!(sdpBasisFull[(mu,1)], K[2][i])

                # if mu == Partition([5,4])
                #     return sym
                # end
            end
        end

        inds = findColSpanningSetV3(aSym)
        if length(inds) > 0

            sdpBasisFull[(mu,-1)] = []

            for i in inds
                push!(sdpBasisFull[(mu,-1)], K[2][i])
            end
        end
    end

    # return blks
    # return test
    return sdpBasisFull#, allReynolds

end

# @time decomposeModule(10)

## Checking if windows basis is correct
# test = load("data/decompose9.jld2")["basis"]

# reynolds = decomposeModule(9)

# # (mu, s) = rand(keys(test))

# # for ((mu,s), bs) in test
# #     R = s == 1 ? reynolds[mu][2] : reynolds[mu][3]
# #     cols =  indexin(test[(mu,s)], reynolds[mu][1][2])
# #     if rank(R[:, cols]) != length(test[(mu,s)])
# #         @error("Incorrect basis!")
# #     end
# # end

# for k in keys(test)
#     @show length(test[k]) == length(reynolds[k])
#     if length(test[k]) != length(reynolds[k])
#         @error("Wrong blocksizes")
#     end
# end