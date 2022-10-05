include("CalcProduct.jl")
include("Decompose.jl")


function blockDiagonalize(basis)
    blockDiag = Dict()
    for ((lambda, sgn), bs) in basis
        m = sum(lambda)
        @show (lambda, sgn, length(bs))
        s = length(bs)
        blockDiag[(lambda, sgn)] = Dict([])
        
        tmpB = zeros(Int, m)
        tmpC = zeros(Int, m)
        for i = 1:s
            for j = i:s
                print("($i, $j)    \n")
                prod = calcProductNoDictThinMat(bs[i],bs[j])
                # @show prod
                invProd = Dict(labelCanonical(Int.(reverse(c)); cAdj = tmpB, best = tmpC) => sgn * coeff for (c,coeff) in prod)
                # @show invProd

                entry = mergewith!(+, prod, invProd)


                # Cannot just merge with -, if keys are different!
                # entry = mergewith!(sgn == 1 ? (+) : (-), prod, invProd)
                filter!(x->x.second != 0, entry )
                # @show entry
                for (cycle, coeff) in entry
                    if !haskey(blockDiag[(lambda, sgn)], cycle)
                        blockDiag[(lambda, sgn)][cycle] = zeros(T, s,s)
                    end
                    blockDiag[(lambda, sgn)][cycle][i,j] = coeff 
                end
            end
        end
    end
    return blockDiag
end

# @time @inbounds test = blockDiagonalize(decomposeModule(8))
# @time CrossingBlockDiagFULL(8)

##

# bs = rand(decomposeModule(10))[2]
# b1 = rand(bs)
# b2 = rand(bs)

# @time @inbounds calcProduct(b1,b2)



# @time unique([labelCanonical(AbstractAlgebra.Generic.matrix_repr(p)) for p in SymmetricGroup(7)])