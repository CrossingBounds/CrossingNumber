# Calculate products of polytabloids and symmetrize them

using AbstractAlgebra, DataStructures

# labelDict = Dict{BitMatrix, Vector{Int8}}()
# labelDictCycle = Dict{Vector{Int8}, Vector{Int8}}()

# function labelCanonical(cycle::Vector{Int8})
#     # return get!(labelDictCycle, cycle) do
#     n = length(cycle)
#     M = BitArray(undef, n, n)
#     M .= 0
#     for i = 1:n
#         M[i, cycle[i]] = 1
#     end
#     return labelCanonical(M, zeros(Int, n), zeros(Int, n), zeros(Int, n))
#     # end
# end

function labelCanonical!(cycle::Vector{T}, cAdj::Vector{T} = zeros(T, length(cycle)), best::Vector{T} = zeros(T, length(cycle))) where {T<:Integer}
    n = T(length(cycle))
    best .= cycle


    function lbl(c::Vector{T})
        t::T = zero(T)
        for i = 1:n
            @inbounds t = c[i]

            for j in eachindex(c)
                @inbounds if c[j] < t
                    @inbounds c[j] += n + T(1) - t
                else
                    @inbounds c[j] += T(1) - t
                end
            end

            isSmaller = true
            for (a, b) in Iterators.zip(Iterators.flatten((i:n, 1:i-1)), Iterators.flatten((1:(n-i+1), n-i+2:n)))
                if @inbounds c[a] < best[b]
                    break
                elseif @inbounds c[a] > best[b]
                    isSmaller = false
                    break
                end
            end

            # @inbounds @views if (c[i:n], c[1:i-1]) < (best[1:(n-i+1)], best[n-i+2:end])
            # @inbounds @views if c[collect(Iterators.flatten((i:n,1:i-1)))] < best[collect(Iterators.flatten((1:(n-i+1),n-i+2:end)))]
            if isSmaller
                @views @inbounds best[1:(n-i+1)] = c[i:n]
                @views @inbounds best[n-i+2:end] = c[1:i-1]
            end
        end
        # return best
        # nothing
        # return minimum(((circshift(cycle,-i) .+ (n-cycle[i+1])) .% n) .+ 1 for i = 0:n-1)
    end

    lbl(cycle)

    # cAdj = zero(cycle)
    for i = 1:n
        @inbounds cAdj[cycle[i]] = i
    end
    lbl(cAdj)

    reverse!(cycle)
    cycle .= n .- cycle .+ 1

    lbl(cycle)
    reverse!(cAdj)
    cAdj .= n .- cAdj .+ 1
    lbl(cAdj)

    # MInv = M[end:-1:1, end:-1:1]
    cycle .= best
    # return Int8.(best)#min(lbl(M), lbl(MInv), lbl(M'), lbl(MInv'))
    # end
end

function labelCanonical(cycle::Vector{Int}; cAdj::Vector{Int} = zeros(Int, length(cycle)), best::Vector{Int} = zeros(Int, length(cycle)))
    n = length(cycle)
    best .= cycle


    function lbl(c::Vector{Int})
        for i = 1:n
            @inbounds t = c[i]

            for j in eachindex(c)
                @inbounds if c[j] < t
                    @inbounds c[j] += n + 1 - t
                else
                    @inbounds c[j] += 1 - t
                end
            end

            isSmaller = true
            for (a, b) in Iterators.zip(Iterators.flatten((i:n, 1:i-1)), Iterators.flatten((1:(n-i+1), n-i+2:n)))
                if @views @inbounds c[a] < best[b]
                    break
                elseif @views @inbounds c[a] > best[b]
                    isSmaller = false
                    break
                end
            end

            # @inbounds @views if (c[i:n], c[1:i-1]) < (best[1:(n-i+1)], best[n-i+2:end])
            # @inbounds @views if c[collect(Iterators.flatten((i:n,1:i-1)))] < best[collect(Iterators.flatten((1:(n-i+1),n-i+2:end)))]
            if isSmaller
                @views @inbounds best[1:(n-i+1)] .= c[i:n]
                @views @inbounds best[n-i+2:end] .= c[1:i-1]
            end
        end
        return best
        # nothing
        # return minimum(((circshift(cycle,-i) .+ (n-cycle[i+1])) .% n) .+ 1 for i = 0:n-1)
    end

    lbl(cycle)

    # cAdj = zero(cycle)
    for i = 1:n
        @inbounds cAdj[cycle[i]] = i
    end
    lbl(cAdj)

    reverse!(cycle)
    cycle .= n .- cycle .+ 1

    lbl(cycle)
    reverse!(cAdj)
    cAdj .= n .- cAdj .+ 1
    lbl(cAdj)

    # MInv = M[end:-1:1, end:-1:1]

    return Int8.(best)#min(lbl(M), lbl(MInv), lbl(M'), lbl(MInv'))
    # end
end

function labelCanonical(M::BitMatrix, cycle::Vector{Int}, cAdj::Vector{Int}, best::Vector{Int})
    # display(M)
    # return get!(labelDict, M) do
    n = Int(size(M, 1))
    # f = (x,t) -> (y = x-t+1; y > 0 ? y : y + n)
    # @inbounds tmp = Int8[0 for i = 1:n]
    # cycle = Int8[n+1 for i = 1:n]
    # cycle .= M * (Int.(1:n))
    for i = 1:n
        @views @inbounds cycle[i] = findfirst(M[:, i])
    end

    labelCanonical(cycle, cAdj, best)
    
end

function orbit(cycle::Vector{Int8})
    n = length(cycle)
    M = BitArray(undef, n, n)
    M .= 0
    for i = 1:n
        M[i, cycle[i]] = 1
    end
    function shifts(M)
        cycle = Int8[findfirst(M[:, i]) for i = 1:n]
        return [((circshift(cycle, -i) .+ (n - cycle[i+1])) .% n) .+ 1 for i = 0:n-1]
    end

    MInv = M[n:-1:1, n:-1:1]

    return unique(vcat(shifts(M), shifts(MInv), shifts(M'), shifts(MInv')))

end

function orbitSize(cycle::Vector{Int8})
    n = length(cycle)
    M = BitArray(undef, n, n)
    M .= 0
    for i = 1:n
        M[i, cycle[i]] = 1
    end
    function shifts(M)
        cycle = Int8[findfirst(M[:, i]) for i = 1:n]
        return [((circshift(cycle, -i) .+ (n - cycle[i+1])) .% n) .+ 1 for i = 0:n-1]
    end

    MInv = M[n:-1:1, n:-1:1]

    return length(unique(vcat(shifts(M), shifts(MInv), shifts(M'), shifts(MInv'))))

end

# with Cm symmetry
# use with @inbounds?
const T = Int64
function calcProduct(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    #TODO: remove fact (since it is always 1)
    # display(t1)
    # display(t2)

    m = max(maximum(t1), maximum(t2))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    @inline function move!(res::Matrix{Int8}, A::Matrix{Int8}, i::Int, j::Int, k::Int, l::Int)
        @inbounds res .= A
        @inbounds res[i, j] -= 1
        @inbounds res[k, l] += 1
        nothing
        # return res
    end

    function pMat(p::Perm{T})
        t = zeros(Int8, m, m)
        for (i, d) in enumerate(p.d)
            t[i, d] = 1
        end
        return t
    end

    function oDet(n::T)
        return Dict(pMat(P) => sign(P) for P in SymmetricGroup(n))
    end

    function mink(As::Dict{Matrix{Int8},T}, Bs::Dict{Matrix{Int8},T})
        t = Dict{Matrix{Int8},T}()#zero(T))
        for (A, ca) in As
            for (B, cb) in Bs
                C = A + B
                t[C] = get!(t, C, 0) + ca * cb
            end
        end
        return t
    end

    res = Dict{Matrix{Int8},T}()#(zero(T))
    res[zeros(Int8, m, m)] = 1
    for k = 1:m
        t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(T(k))
            map!(x -> x * factorial(k), values(oD))
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    B2 = zeros(Int8, m, m)

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                res2 = Dict{Matrix{Int8},T}()#(zero(T))
                for (B, c) in res
                    for nz in findall(!iszero, B[:, j])
                        move!(B2, B, nz, j, nz, s)
                        kB2 = copy(B2)
                        @inbounds res2[kB2] = get!(res2, kB2, 0) + c * B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                end
                res = res2
            end
            # @show length(res)
        end
    end

    # @show length(res)
    res2 = Dict{Matrix{Int8},T}()#(zero(T))
    for (B, c) in res
        # display(B)
        B3 = reshape(maximum([vec(circshift(B, (0, i))) for i = 0:size(B, 1)-1]), size(B))
        res2[B3] = get!(res2, B3, 0) + c
        # res2[B3] += c
    end
    filter!(x -> x[2] != 0, res2)
    res = res2
    # @show length(res)


    # @info("Loop 2")

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                res2 = Dict{Matrix{Int8},T}()#(zero(T))
                for (B, c) in res
                    for nz in findall(!iszero, B[j, :])
                        move!(B2, B, j, nz, s, nz)
                        kB2 = copy(B2)
                        @inbounds res2[kB2] = get!(res2, kB2, 0) + c * B[j, nz]

                        # res2[copy(B2)] += c*B[j,nz]
                    end
                end
                res = res2
            end
            # @show length(res)
        end
    end

    # @show fact
    map!(x -> Int(x / fact), values(res))

    # tmpA = zeros(Int8, m)
    # tmpB = zeros(Int8, m)

    labelledRes = Dict{Vector{Int8},T}()
    for (M, c) in res
        labelledC = labelCanonical(Bool.(M))#, tmpA, tmpB)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

function calcProductNoDict(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    #TODO: remove fact (since it is always 1)
    # display(t1)
    # display(t2)

    m = max(maximum(t1), maximum(t2))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    @inline function move!(res::Matrix{Int8}, A::Vector{Int8}, i::Int, j::Int, k::Int, l::Int)
        res .= reshape(A, m, m)
        @inbounds res[i, j] -= 1
        @inbounds res[k, l] += 1
        nothing
        # return res
    end

    function pMat(p::Perm{T})
        t = zeros(Int8, m, m)
        for (i, d) in enumerate(p.d)
            @inbounds t[i, d] = 1
        end
        return t
    end

    function oDet(n::T)
        return [(vec(pMat(P)), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Vector{Int8},T}}, Bs::Vector{Tuple{Vector{Int8},T}})
        t = Vector{Tuple{Vector{Int8},T}}()#zero(T))
        for (A, ca) in As
            for (B, cb) in Bs
                C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                push!(t, (C, ca * cb))
            end
        end
        return t
    end

    # function reduceVec(v)
    #     # return v
    #     preL = length(v)
    #     @inbounds @views sort!(v, by = t -> t[1])
    #     tmp = Vector{Tuple{Vector{Int8},T}}()
    #     last = nothing
    #     s = 0
    #     for (M, c) in v
    #         if M == last
    #             s += c
    #         else
    #             if last !== nothing && s != 0
    #                 push!(tmp, (last, s))
    #             end
    #             last = M
    #             s = c
    #         end
    #     end
    #     if last !== nothing && s != 0
    #         push!(tmp, (last, s))
    #     end
    #     # ratio = length(tmp)/preL
    #     # ratio < 1 && @info("Reduced $(ratio)")
    #     return tmp
    # end

    # in place
    function reduceVec!(v)
        # return v
        # @show v
        preL = length(v)
        @inbounds @views sort!(v, by = t -> t[1])
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        return v
    end

    res = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, (vec(zeros(Int8, m, m)), 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(T(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # res = reduceVec(res)
    reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    B2 = zeros(Int8, m, m)

    tmp = zeros(Int8, m, m)

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                for (B, c) in res
                    tmp = reshape(B, m, m)
                    @inbounds @views for nz in findall(!iszero, tmp[:, j])
                        @inbounds move!(B2, B, nz, j, nz, s)
                        kB2 = copy(vec(B2))
                        @inbounds push!(res2, (kB2, c * tmp[nz, j]))
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                end
                res = reduceVec!(res2)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    for (B, c) in res
        # display(B)
        @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        push!(res2, (B3, c))
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    res = reduceVec!(res2)
    # @show length(res)


    # @info("Loop 2")

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                for (B, c) in res
                    @inbounds tmp = reshape(B, m, m)
                    @views @inbounds for nz in findall(!iszero, tmp[j, :])
                        @inbounds move!(B2, B, j, nz, s, nz)
                        kB2 = copy(vec(B2))
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        @inbounds push!(res2, (kB2, c * tmp[j, nz]))

                        # res2[copy(B2)] += c*B[j,nz]
                    end
                end
                res = reduceVec!(res2)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    labelledRes = Dict{Vector{Int8},T}()
    for (M, c) in res
        @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
        labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

using ProgressMeter

function calcProductNoDict2(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    m = max(maximum(t1), maximum(t2))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    function move!(res::Matrix{Int8}, A::Vector{Int8}, i::Int, j::Int, k::Int, l::Int)
        @views res .= reshape(A, m, m)
        @inbounds res[i, j] -= 1
        @inbounds res[k, l] += 1
        nothing
    end

    function pMat(p::Perm{T})
        t = zeros(Int8, m, m)
        for (i, d) in enumerate(p.d)
            @inbounds t[i, d] = 1
        end
        return t
    end

    function oDet(n::T)
        return [(vec(pMat(P)), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Vector{Int8},T}}, Bs::Vector{Tuple{Vector{Int8},T}})
        t = Vector{Tuple{Vector{Int8},T}}()#zero(T))
        for (A, ca) in As
            for (B, cb) in Bs
                C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                push!(t, (C, ca * cb))
            end
        end
        return t
    end

    # in place
    function reduceVec!(v)
        # return v
        # @show v
        preL = length(v)
        @inbounds @views sort!(v, by = t -> t[1])
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        return v
    end

    res = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, (vec(zeros(Int8, m, m)), 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(T(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # res = reduceVec(res)
    reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    B2 = zeros(Int8, m, m)

    tmp = zeros(Int8, m, m)


    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @inbounds tmp .= reshape(B, m, m)
                    isFirst = true
                    @inbounds @views for nz in findall(!iszero, tmp[:, j])
                        @inbounds move!(B2, B, nz, j, nz, s)
                        @views kB2 = copy(vec(B2))
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[nz, j])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[nz, j]))
                        end
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                    isFirst && @error("Should not be first")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp2 = zero(tmp)
    best = zero(tmp)

    for (i, (B, c)) in enumerate(res)
        tmp .= reshape(B, m, m)
        best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, (0, j))
            if @views vec(tmp2) < vec(best)
                @inbounds best .= tmp2
            end
        end

        @inbounds @views B .= vec(best)

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                @showprogress for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @views @inbounds tmp .= reshape(B, m, m)

                    isFirst = true
                    @views @inbounds for nz in findall(!iszero, tmp[j, :])
                        @inbounds move!(B2, B, j, nz, s, nz)
                        @views kB2 = copy(vec(B2))
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[j, nz])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[j, nz]))
                        end
                        # res2[copy(B2)] += c*B[j,nz]
                    end
                    isFirst && @error("Should not be first")
                end
                # res = reduceVec!(res2)
                reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    tmpA = zeros(Int, m)
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)

    labelledRes = Dict{Vector{Int8},T}()
    for (M, c) in res
        labelledC = labelCanonical(Bool.(reshape(M, m, m)), tmpA, tmpB, tmpC)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    # labelledRes = Dict{Vector{Int8},T}()
    # for (M, c) in res
    #     @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
    #     labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    # end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

function calcProductNoDictVec(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    m = max(maximum(t1), maximum(t2))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    function move!(res::Vector{Int8}, A::Vector{Int8}, i::Int, j::Int, k::Int, l::Int)
        @views res .= A#reshape(A, m, m)
        @inbounds res[i+(j-1)*m] -= 1
        @inbounds res[k+(l-1)*m] += 1
        nothing
    end

    function pMat(p::Perm{Int})
        # t = zeros(Int8, m, m)
        t = zeros(Int8, m^2)
        for (i, d) in enumerate(p.d)
            @inbounds t[i+(d-1)*m] = 1
        end
        return t
    end

    function oDet(n::Int)
        return [(pMat(P), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Vector{Int8},Int}}, Bs::Vector{Tuple{Vector{Int8},Int}})
        t = Vector{Tuple{Vector{Int8},Int}}()#zero(T))
        for (A, ca) in As
            for (B, cb) in Bs
                C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                push!(t, (C, ca * cb))
            end
        end
        return t
    end

    # in place
    function reduceVec!(v)
        # @show length(v)
        # display(varinfo(all = true, sortby = :size))
        # @show Base.summarysize(v)
        # @show Base.summarysize(v)/length(v)
        preL = length(v)


        # @show v

        # @show issorted(v)

        # function HoarePartition(A, lo, hi)
        #     @inbounds pivot = A[floor(Int, (hi+lo)/2)][1]
        #     i = lo - 1
        #     j = hi + 1
        #     while true
        #         i += 1
        #         @inbounds while A[i][1] < pivot 
        #             i += 1
        #         end
        #         j -= 1
        #         @inbounds while A[j][1] > pivot
        #             j -= 1
        #         end
        #         if i >= j
        #             return j
        #         end

        #         @inbounds swap = A[i]
        #         @inbounds A[i] = A[j]
        #         @inbounds A[j] = swap
        #     end
        # end

        # function HoareQuickSort(A, lo = 1, hi = length(A))
        #     if lo >= 0 && hi >= 0 && lo < hi 
        #         p = HoarePartition(A, lo, hi)
        #         HoareQuickSort(A, lo, p)
        #         HoareQuickSort(A, p+1, hi)
        #     end
        # end

        # HoareQuickSort(v)

        @inbounds @views sort!(v, by = t -> t[1])

        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        # @show v
        # @show (preL, length(v))
        return v
    end

    res = Vector{Tuple{Vector{Int8},Int}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, (zeros(Int8, m^2), 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(Int(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # res = reduceVec(res)
    # reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    # B2 = zeros(Int8, m, m)
    B2 = zeros(Int8, m^2)

    # tmp = zeros(Int8, m, m)
    tmp = zeros(Int8, m^2)


    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @inbounds tmp .= B#reshape(B, m, m)
                    isFirst = true
                    # @inbounds @views for nz in findall(!iszero, tmp[:, j])

                    @inbounds @views for nz in findall(!iszero, tmp[(1+(j-1)*m):(m+(j-1)*m)])
                        @inbounds move!(B2, B, nz, j, nz, s)
                        # @views kB2 = copy(vec(B2))
                        kB2 = copy(B2)
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[nz+(j-1)*m])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[nz+(j-1)*m]))
                        end
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                    isFirst && @error("Should not be first A")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                # reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp2 = zero(tmp)
    best = zero(tmp)

    for (i, (B, c)) in enumerate(res)
        @inbounds tmp .= B#reshape(B, m, m)
        @inbounds best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, j * m)
            if @views vec(tmp2) < vec(best)
                @inbounds best .= tmp2
            end
        end

        @inbounds B .= best

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                @showprogress 1 "Innermost loop..." for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @views @inbounds tmp .= B#reshape(B, m, m)

                    isFirst = true
                    # @views @inbounds for nz in findall(!iszero, tmp[j, :])
                    @views @inbounds for nz in findall(!iszero, tmp[(j):m:(j+(m-1)*m)])
                        @inbounds move!(B2, B, j, nz, s, nz)
                        kB2 = copy(B2)
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[j+(nz-1)*m])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[j+(nz-1)*m]))
                        end
                        # res2[copy(B2)] += c*B[j,nz]
                    end
                    isFirst && @error("Should not be first B")
                end
                # res = reduceVec!(res2)
                # @info("Reducing...")
                # reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    tmpA = zeros(Int, m)
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)

    labelledRes = Dict{Vector{Int8},Int}()
    @showprogress 1 "labeling canonically..." for (M, c) in res
        labelledC = labelCanonical(Bool.(reshape(M, m, m)), tmpA, tmpB, tmpC)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    # labelledRes = Dict{Vector{Int8},T}()
    # for (M, c) in res
    #     @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
    #     labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    # end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

using SparseArrays

## different way to represent overlap matrices: sparse vectorized matrix (SLOW :c ). 
function calcProductNoDictVecSparse(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    m = max(maximum(t1), maximum(t2))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    function move!(res::SparseVector{Int8,Int}, A::SparseVector{Int8,Int}, i::Int, j::Int, k::Int, l::Int)
        @views res .= A#reshape(A, m, m)
        @inbounds res[i+(j-1)*m] -= 1
        @inbounds res[k+(l-1)*m] += 1
        dropzeros!(res)
        nothing
    end

    function pMat(p::Perm{Int})
        # t = zeros(Int8, m, m)
        t = spzeros(Int8, m^2)
        for (i, d) in enumerate(p.d)
            @inbounds t[i+(d-1)*m] = 1
        end
        return t
    end

    function oDet(n::Int)
        return [(pMat(P), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{SparseVector{Int8,Int},Int}}, Bs::Vector{Tuple{SparseVector{Int8,Int},Int}})
        t = Vector{Tuple{SparseVector{Int8,Int},Int}}()#zero(T))
        for (A, ca) in As
            for (B, cb) in Bs
                C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                push!(t, (C, ca * cb))
            end
        end
        return t
    end

    # in place
    function reduceVec!(v)
        preL = length(v)
        function HoarePartition(A, lo, hi)
            @inbounds pivot = A[floor(Int, (hi + lo) / 2)][1]
            i = lo - 1
            j = hi + 1
            while true
                i += 1
                @inbounds while A[i][1] < pivot
                    i += 1
                end
                j -= 1
                @inbounds while A[j][1] > pivot
                    j -= 1
                end
                if i >= j
                    return j
                end

                @inbounds swap = A[i]
                @inbounds A[i] = A[j]
                @inbounds A[j] = swap
            end
        end

        function HoareQuickSort(A, lo = 1, hi = length(A))
            if lo >= 0 && hi >= 0 && lo < hi
                p = HoarePartition(A, lo, hi)
                HoareQuickSort(A, lo, p)
                HoareQuickSort(A, p + 1, hi)
            end
        end

        HoareQuickSort(v)
        # @inbounds @views sort!(v, by = t -> t[1])

        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        # @show v
        return v
    end

    res = Vector{Tuple{SparseVector{Int8,Int},Int}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, (spzeros(Int8, m^2), 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(Int(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # res = reduceVec(res)
    reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    # B2 = zeros(Int8, m, m)
    B2 = spzeros(Int8, m^2)

    # tmp = zeros(Int8, m, m)
    tmp = spzeros(Int8, m^2)


    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @inbounds tmp .= B#reshape(B, m, m)
                    isFirst = true
                    # @inbounds @views for nz in findall(!iszero, tmp[:, j])

                    @inbounds @views for nz in findall(!iszero, tmp[(1+(j-1)*m):(m+(j-1)*m)])
                        @inbounds move!(B2, B, nz, j, nz, s)
                        # @views kB2 = copy(vec(B2))
                        kB2 = copy(B2)
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[nz+(j-1)*m])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[nz+(j-1)*m]))
                        end
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                    isFirst && @error("Should not be first A")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp2 = zero(tmp)
    best = zero(tmp)

    for (i, (B, c)) in enumerate(res)
        @inbounds tmp .= B#reshape(B, m, m)
        @inbounds best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, j * m)
            if @views vec(tmp2) < vec(best)
                @inbounds best .= tmp2
            end
        end

        @inbounds B .= best

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-1:-1:1
        # @show j
        # First loop depends on order, second does not matter
        for s = m:-1:j+1
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                @showprogress 1 "Innermost loop..." for i in 1:curL
                    @inbounds (B, c) = res[i]
                    @views @inbounds tmp .= B#reshape(B, m, m)

                    isFirst = true
                    # @views @inbounds for nz in findall(!iszero, tmp[j, :])
                    @views @inbounds for nz in findall(!iszero, tmp[(j):m:(j+(m-1)*m)])
                        @inbounds move!(B2, B, j, nz, s, nz)
                        kB2 = copy(B2)
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        if isFirst
                            @inbounds res[i] = (kB2, c * tmp[j+(nz-1)*m])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * tmp[j+(nz-1)*m]))
                        end
                        # res2[copy(B2)] += c*B[j,nz]
                    end
                    isFirst && @error("Should not be first B")
                end
                # res = reduceVec!(res2)
                # @info("Reducing...")
                reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    tmpA = zeros(Int, m)
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)

    labelledRes = Dict{Vector{Int8},Int}()
    @showprogress 1 "labeling canonically..." for (M, c) in res
        labelledC = labelCanonical(Bool.(Matrix(reshape(M, m, m))), tmpA, tmpB, tmpC)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    # labelledRes = Dict{Vector{Int8},T}()
    # for (M, c) in res
    #     @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
    #     labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    # end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

## different way to represent overlap matrices: Vector{Vector{Int8}}. 
function calcProductNoDictVecVec(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    m = Int8(max(maximum(t1), maximum(t2)))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    # move entry from (i,j) -> (k,l)
    function move!(res::Vector{Vector{Int8}}, A::Vector{Vector{Int8}}, i::Int8, j::Int8, k::Int8, l::Int8)
        # @show length(vcat(A...))
        @inbounds res .= A#reshape(A, m, m)
        # @show A
        # @show (i, j, k, l)
        res[j] = copy(res[j])
        res[l] = copy(res[l])
        deleteat!(res[j], findfirst(x -> x == i, res[j]))
        push!(res[l], k)
        # @show length(vcat(res...))

        # if length(vcat(A...)) != length(vcat(res...))
            # @error("Something wrogn")
        # end
        # sort!()
        # @inbounds res[i+(j-1)*m] -= 1
        # @inbounds res[k+(l-1)*m] += 1
        # dropzeros!(res)
        nothing
    end

    function pMat(p::Perm{Int8})
        # t = zeros(Int8, m, m)
        # t = spzeros(Int8, m^2)
        # for (i, d) in enumerate(p.d)
        #     @inbounds t[i+(d-1)*m] = 1
        # end
        return vcat([[d] for d in p.d], [Int8[] for i = length(p.d)+1:m])
    end

    function oDet(n::Int8)
        return [(pMat(P), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Vector{Vector{Int8}},Int}}, Bs::Vector{Tuple{Vector{Vector{Int8}},Int}})
        t = Vector{Tuple{Vector{Vector{Int8}},Int}}()#zero(T))
        # @show As
        # @show Bs
        for (A, ca) in As
            for (B, cb) in Bs
                # C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                C = [sort!(vcat(a, b)) for (a, b) in zip(A, B)]
                push!(t, (C, ca * cb))
            end
        end
        # @show t
        return t
    end

    # in place
    function reduceVec!(v)
        preL = length(v)
        
        @inbounds @views sort!(v, by = t -> t[1])

        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        # @show v
        return v
    end

    res = Vector{Tuple{Vector{Vector{Int8}},Int}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, ([Int8[] for i = 1:m], 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(Int8(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            # @show k
            # @show oD
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # @show res

    # res = reduceVec(res)
    reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    # B2 = zeros(Int8, m, m)
    # B2 = spzeros(Int8, m^2)
    B2 = [Int8[] for i = 1:m]

    # tmp = zeros(Int8, m, m)
    # tmp = spzeros(Int8, m^2)
    # tmp = [Int8[] for i = 1:m]


    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]

                    # @show B
                    # @show typeof(B)
                    # @inbounds tmp .= B#reshape(B, m, m)
                    isFirst = true
                    # @inbounds @views for nz in findall(!iszero, tmp[:, j])

                    @inbounds @views for nz in unique(B[j])
                        # @show nz
                        # @show typeof(nz)
                        # kB2 = [Int8[] for i = 1:m]
                        # @inbounds move!(kB2, B, nz, j, nz, s)
                        # @views kB2 = copy(vec(B2))
                        # kB2 = copy(B2)

                        kB2 = copy(B)
                        @inbounds kB2[j] = copy(B[j])
                        @inbounds kB2[s] = copy(B[s])
                        @inbounds deleteat!(kB2[j], findfirst(x -> x == nz, kB2[j]))
                        @inbounds push!(kB2[s], nz)

                        if isFirst
                            @inbounds res[i] = (kB2, c * count(x -> x == nz, B[j]))
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * count(x -> x == nz, B[j])))
                        end
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                    isFirst && @error("Should not be first A")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp = [Int8[] for i = 1:m]
    tmp2 = [Int8[] for i = 1:m]
    best = [Int8[] for i = 1:m]

    for (i, (B, c)) in enumerate(res)
        @inbounds tmp .= B#reshape(B, m, m)
        @inbounds best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, j)
            if @views tmp2 < best
                @inbounds best .= tmp2
            end
        end

        # transpose
        for j = 1:m 
            B[j] = Int8[]
        end
        for j = Int8(1):m 
            for d in best[j]
                push!(B[d], j)
            end
        end

        # @inbounds B .= best

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                # @showprogress 1 "Innermost loop..." 
                for i in 1:curL
                    @inbounds (B, c) = res[i]
                    # @views @inbounds tmp .= B#reshape(B, m, m)

                    isFirst = true
                    # @views @inbounds for nz in findall(!iszero, tmp[j, :])
                    # @views @inbounds for nz in findall(!iszero, tmp[(j):m:(j+(m-1)*m)])
                    for nz in unique(B[j])
                        # kB2 = [Int8[] for i = 1:m]
                        # @inbounds move!(kB2, B, nz, j, nz, s)
                        
                        kB2 = copy(B)
                        @inbounds kB2[j] = copy(B[j])
                        @inbounds kB2[s] = copy(B[s])
                        @inbounds deleteat!(kB2[j], findfirst(x -> x == nz, kB2[j]))
                        @inbounds push!(kB2[s], nz)

                        # kB2 = copy(B2)
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        if isFirst
                            @inbounds res[i] = (kB2, c * count(x -> x == nz, B[j]))
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * count(x -> x == nz, B[j])))
                        end
                        # res2[copy(B2)] += c*B[j,nz]
                    end
                    isFirst && @error("Should not be first B")
                end
                # res = reduceVec!(res2)
                # @info("Reducing...")
                reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    cycle = zeros(Int, m)
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)

    # tmpM = BitMatrix(zeros(Bool, m, m))
    # M = BitArray(undef, m, m)
    
    labelledRes = Dict{Vector{Int8},Int}()
    # @showprogress 1 "labeling canonically..." 
    for (B, c) in res
        # M .= 0
        for i = 1:m
            # M[i, B[i][1]] = 1
            cycle[i] = B[i][1]
        end
        
        labelledC = labelCanonical(cycle, tmpB, tmpC)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    # labelledRes = Dict{Vector{Int8},T}()
    # for (M, c) in res
    #     @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
    #     labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    # end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end

function calcProductNoDictVecVec2(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

    m = Int8(max(maximum(t1), maximum(t2)))


    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    # move entry from (i,j) -> (k,l)
    function move!(res::Vector{Vector{Int8}}, A::Vector{Vector{Int8}}, i::Int8, j::Int8, k::Int8, l::Int8)
        @inbounds res .= A
        # @show A
        # @show (i, j, k, l)
        res[j] = copy(res[j])
        res[l] = copy(res[l])
        # deleteat!(res[j], findfirst(x -> x == i, res[j]))
        res[j][i] -= 1
        # push!(res[l], k)
        res[l][k] += 1
        # @show length(vcat(res...))

        # if length(vcat(A...)) != length(vcat(res...))
            # @error("Something wrogn")
        # end
        # sort!()
        # @inbounds res[i+(j-1)*m] -= 1
        # @inbounds res[k+(l-1)*m] += 1
        # dropzeros!(res)
        nothing
    end

    function pMat(p::Perm{Int8})
        # t = zeros(Int8, m, m)
        # t = spzeros(Int8, m^2)
        # for (i, d) in enumerate(p.d)
        #     @inbounds t[i+(d-1)*m] = 1
        # end
        return vcat([[j == d ? Int8(1) : Int8(0) for j = 1:m] for d in p.d], [Int8[0 for j = 1:m] for i = length(p.d)+1:m])
    end

    function oDet(n::Int8)
        return [(pMat(P), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Vector{Vector{Int8}},Int}}, Bs::Vector{Tuple{Vector{Vector{Int8}},Int}})
        t = Vector{Tuple{Vector{Vector{Int8}},Int}}()#zero(T))
        # @show As
        # @show Bs
        for (A, ca) in As
            for (B, cb) in Bs
                # C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                # C = [sort!(vcat(a, b)) for (a, b) in zip(A, B)]
                C = [a + b for (a, b) in zip(A, B)]
                
                push!(t, (C, ca * cb))
            end
        end
        # @show t
        return t
    end

    # in place
    function reduceVec!(v)
        preL = length(v)
        
        @inbounds @views sort!(v, by = t -> t[1])

        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        # @show v
        return v
    end

    res = Vector{Tuple{Vector{Vector{Int8}},Int}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, ([Int8[0 for j = 1:m] for i = 1:m], 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(Int8(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            # @show k
            # @show oD
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # @show res

    # res = reduceVec(res)
    reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    # B2 = zeros(Int8, m, m)
    # B2 = spzeros(Int8, m^2)
    B2 = [Int8[0 for j = 1:m] for i = 1:m]

    # tmp = zeros(Int8, m, m)
    # tmp = spzeros(Int8, m^2)
    # tmp = [Int8[] for i = 1:m]


    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]

                    # @show B
                    # @show typeof(B)
                    # @inbounds tmp .= B#reshape(B, m, m)
                    isFirst = true
                    # @inbounds @views for nz in findall(!iszero, tmp[:, j])

                    @inbounds @views for nz::Int8 in findall(B[j] .!= 0)
                        # @show nz
                        # @show typeof(nz)
                        @inbounds move!(B2, B, nz, j, nz, s)
                        # @views kB2 = copy(vec(B2))
                        kB2 = copy(B2)
                        if isFirst
                            @inbounds res[i] = (kB2, c * B[j][nz])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * B[j][nz]))
                        end
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                        # res2[copy(B2)] += c*B[nz, j]
                    end
                    isFirst && @error("Should not be first A")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res

    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp = [Int8[0 for j = 1:m] for i = 1:m]
    tmp2 = [Int8[0 for j = 1:m] for i = 1:m]
    best = [Int8[0 for j = 1:m] for i = 1:m]

    for (i, (B, c)) in enumerate(res)
        @inbounds tmp .= B#reshape(B, m, m)
        @inbounds best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, j)
            if @views tmp2 < best
                @inbounds best .= tmp2
            end
        end

        # transpose
        for j = 1:m 
            B[j] = Int8[0 for o = 1:m]
        end
        for j = Int8(1):m 
            for d = 1:m
            # for d in best[j]
                # push!(B[d], j)
                B[d][j] = best[j][d]
            end
        end

        # @inbounds B .= best

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                @showprogress 1 "Innermost loop..." for i in 1:curL
                    @inbounds (B, c) = res[i]
                    # @views @inbounds tmp .= B#reshape(B, m, m)

                    isFirst = true
                    # @views @inbounds for nz in findall(!iszero, tmp[j, :])
                    # @views @inbounds for nz in findall(!iszero, tmp[(j):m:(j+(m-1)*m)])
                    for nz::Int8 in findall(B[j] .!= 0)
                        @inbounds move!(B2, B, nz, j, nz, s)
                        kB2 = copy(B2)
                        # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[j,nz]
                        if isFirst
                            @inbounds res[i] = (kB2, c * B[j][nz])
                            isFirst = false
                        else
                            @inbounds push!(res, (kB2, c * B[j][nz]))
                        end
                        # res2[copy(B2)] += c*B[j,nz]
                    end
                    isFirst && @error("Should not be first B")
                end
                # res = reduceVec!(res2)
                # @info("Reducing...")
                reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))

    tmpA = zeros(Int, m)
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)

    # tmpM = BitMatrix(zeros(Bool, m, m))
    M = BitArray(undef, m, m)
    
    labelledRes = Dict{Vector{Int8},Int}()
    @showprogress 1 "labeling canonically..." for (B, c) in res
        # M .= 0
        # for i = 1:m
        #     M[i, B[i][1]] = 1
        # end
        for i = 1:m 
            for j = 1:m 
                M[i,j] = B[i][j]
            end
        end
        # M .= hcat(B...)
        
        labelledC = labelCanonical(M, tmpA, tmpB, tmpC)
        labelledRes[labelledC] = get(labelledRes, labelledC, 0) + c
    end

    # labelledRes = Dict{Vector{Int8},T}()
    # for (M, c) in res
    #     @inbounds labelledC = labelCanonical(Bool.(reshape(M, m, m)))
    #     labelledRes[labelledC] = get!(labelledRes, labelledC, 0) + c / fact
    # end

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    return labelledRes
end


## different way to represent overlap matrices: 2 x m matrix 
# using StaticArrays
function calcProductNoDictThinMat(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64) 

    m = Int8(max(maximum(t1), maximum(t2)))

    lambda = t1.part
    lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

    function pMat(p::Perm{Int8})
        # return vcat([[d] for d in p.d], [Int8[] for i = length(p.d)+1:m])
        P = zeros(Int8, 2, length(p.d))
        for i = 1:length(p.d)
            P[1,i] = i 
            P[2,i] = p.d[i]
        end
        # display(p)
        # display(P)
        return P
    end

    function oDet(n::Int8)
        return [(pMat(P), sign(P)) for P in SymmetricGroup(n)]
    end

    function mink(As::Vector{Tuple{Matrix{Int8},Int}}, Bs::Vector{Tuple{Matrix{Int8},Int}})
        t = Vector{Tuple{Matrix{Int8},Int}}()#zero(T))
        # @show As
        # @show Bs
        for (A, ca) in As
            for (B, cb) in Bs
                # C = A + B
                # t[C] = get!(t, C, 0) + ca*cb
                # C = [hcat(a, b) for (a, b) in zip(A, B)]
                C = hcat(A, B)
                push!(t, (C, ca * cb))
            end
        end
        # @show t
        return t
    end

    # in place
    function reduceVec!(v)
        preL = length(v)

        # u1 = length(unique([x[1] for x in v]))
        # for i = 1:preL
        #     @inbounds (B, c) = v[i]
        #     # @show B
        #     # display(B)
        #     # @show c
        #     old = copy(B)
        #     @inbounds v[i][1] .= sortslices(B, dims = 2)#,c)
        #     # if old != v[i][1]
        #     #     @info("Different")
        #     #     display(old)
        #     #     display(v[i][1])
        #     # end
        #     # display(v[i][1])
        # end
        
        # u2 = length(unique([x[1] for x in v]))

        # @show (u1, u2)
        


        @inbounds @views sort!(v, by = t -> t[1][:])

        deleteInd = BitVector(undef, preL)
        deleteInd .= true
        # tmp = Vector{Tuple{Vector{Int8},T}}()
        last = nothing
        s = 0
        firstOccurance = 1
        for i in eachindex(v)
            @inbounds (M, c) = v[i]
            if M == last
                s += c
            else
                if last !== nothing
                    # push!(tmp, (last, s))
                    if s == 0
                        @inbounds deleteInd[firstOccurance] = true
                    else
                        @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
                    end
                end
                last = M
                s = c
                firstOccurance = i
                deleteInd[i] = false
            end
        end
        if last !== nothing
            # push!(tmp, (last, s))
            if s == 0
                @inbounds deleteInd[firstOccurance] = true
            else
                @inbounds v[firstOccurance] = (v[firstOccurance][1], s)
            end
        end

        deleteat!(v, deleteInd)
        # ratio = length(tmp)/preL
        # ratio < 1 && @info("Reduced $(ratio)")
        # @show v

        # @show (preL, length(v))
        return v
    end

    res = Vector{Tuple{Matrix{Int8},Int}}()#(zero(T))
    # res[zeros(Int8,m,m)] = 1
    push!(res, (zeros(Int8, 2, 0), 1))
    for k = 1:m
        @inbounds t = lambda[k] - lambda[k+1]
        if t > 0
            oD = oDet(Int8(k))
            @inbounds oD = map(x -> (x[1], x[2] * factorial(k)), oD)
            # @show k
            # @show oD
            for i = 1:t
                res = mink(res, oD)
            end
        end
    end

    # @show res

    # res = reduceVec(res)
    # reduceVec!(res)
    # @show res

    # res = Dict(collect(res)[pLi])
    # display(res)

    limit = false
    r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
    u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


    # factor out the size of the column stabilizer
    fact = prod([factorial(t) for t in conj(t1.part).part])
    # fact = 1 

    # @info("Loop 1")

    B2 = zeros(Int8, 2, m)
    # B2 = spzeros(Int8, m^2)
    # B2 = [Int8[] for i = 1:m]

    # tmp = zeros(Int8, m, m)
    # tmp = spzeros(Int8, m^2)
    # tmp = [Int8[] for i = 1:m]

    covered = BitVector([false for i = 1:m])

    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            usj = u(s, j)

            fact *= factorial(usj)

            for k = 1:usj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                for i in 1:curL
                    @inbounds (B, c) = res[i]

                    # @show (j, s)
                    # display(B)

                    # @show B
                    # @show typeof(B)
                    # @inbounds tmp .= B#reshape(B, m, m)
                    isFirst = true
                    # @inbounds @views for nz in findall(!iszero, tmp[:, j])
                    # @inbounds @views for nz in unique(B[j])
                    covered .= false
                    for col = 1:m
                        nz = B[1,col]
                        @inbounds if B[2, col] == j && !covered[nz]
                            covered[nz] = true

                            kB2 = copy(B)
                            # @inbounds kB2[j] = copy(B[j])
                            # @inbounds kB2[s] = copy(B[s])
                            # @inbounds deleteat!(kB2[j], findfirst(x -> x == nz, kB2[j]))
                            # @inbounds push!(kB2[s], nz)
                            kB2[2, col] = s
                            num = 1
                            for index in col+1:m
                                @inbounds if B[2, index] == j && B[1, index] == nz
                                    num += 1
                                end
                            end

                            if isFirst
                                @inbounds res[i] = (kB2, c * num)
                                isFirst = false
                            else
                                @inbounds push!(res, (kB2, c * num))
                            end
                            # @inbounds res2[kB2] = get!(res2, kB2, 0) + c*B[nz, j]

                            # res2[copy(B2)] += c*B[nz, j]
                        end
                    end
                    isFirst && @error("Should not be first A")
                    # if isFirst
                    #     res[i] = (B, 0)
                    # end
                end
                # res = reduceVec!(res2)
                # reduceVec!(res)
            end
            # @show length(res)
        end
    end


    # @show res
    # reduceVec!(res)
    # @show length(res)
    # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
    tmp = [Int8[] for i = 1:m]
    tmp2 = [Int8[] for i = 1:m]
    best = [Int8[] for i = 1:m]
    B = [Int8[] for i = 1:m]

    for (i, (D, c)) in enumerate(res)
        # @info("Label")
        # display(D)
        empty!.(B)
        for j = 1:m
            @inbounds push!(B[D[2,j]], D[1,j])
        end

        # if !all(length.(B) .== 1)
        #     display(D)
        #     display(B)
        #     @error("Error somewhere")
        #     return nothing
        # end

        @inbounds tmp .= B#reshape(B, m, m)
        @inbounds best .= tmp
        for j = 1:m-1
            @inbounds circshift!(tmp2, tmp, j)
            if @views tmp2 < best
                @inbounds best .= tmp2
            end
        end

        # transpose
        # for j = 1:m 
        #     B[j] = Int8[]
        # end
        # for j = Int8(1):m 
        #     for d in best[j]
        #         push!(B[d], j)
        #     end
        # end

        col = 1
        for j = 1:m 
            for d in best[j]
                @inbounds D[1, col] = j
                @inbounds D[2, col] = d
                col += 1 
            end
        end 


        # @inbounds B .= best

        # display(B)
        # @inbounds B3 = maximum([vec(circshift(reshape(B, m, m), (0, i))) for i = 0:m-1])
        # res2[B3] = get!(res2, B3, 0) + c
        # push!(res2, (B3, c))
        # res[i] = (B3, c)
        # res2[B3] += c
        # @info("Labelled D")
        # display(D)
    end
    # filter!(x->x[2] != 0, res2)
    # res = reduceVec!(res2)
    reduceVec!(res)
    # @show length(res)


    # @info("Loop 2")

    for j = m-Int8(1):Int8(-1):Int8(1)
        # @show j
        # First loop depends on order, second does not matter
        for s = m:Int8(-1):j+Int8(1)
            # @show s
            # print("($s,$j)\r")

            rsj = r(s, j)

            fact *= factorial(rsj)


            # res = unique(res)

            for k = 1:rsj
                # res2 = Vector{Tuple{Vector{Int8},T}}()#(zero(T))
                curL = length(res)
                @showprogress 1 "Innermost loop..." for i in 1:curL
                    @inbounds (B, c) = res[i]
                    # @views @inbounds tmp .= B#reshape(B, m, m)

                    isFirst = true
                    # @views @inbounds for nz in findall(!iszero, tmp[j, :])
                    # @views @inbounds for nz in findall(!iszero, tmp[(j):m:(j+(m-1)*m)])
                    covered .= false
                    for col = 1:m
                        nz = B[1,col]
                        @inbounds if B[2, col] == j && !covered[nz]
                            covered[nz] = true

                            kB2 = copy(B)
                            kB2[2, col] = s

                            num = 1
                            for index in col+1:m
                                @inbounds if B[2, index] == j  && B[1, index] == nz
                                    num += 1
                                end
                            end

                            if isFirst
                                @inbounds res[i] = (kB2, c * num)
                                isFirst = false
                            else
                                @inbounds push!(res, (kB2, c * num))
                            end
                        end
                    end
                    
                    isFirst && @error("Should not be first B")
                end
                # res = reduceVec!(res2)
                # @info("Reducing...")
                # reduceVec!(res)
            end
            # @show length(res)
        end
    end

    # @show fact
    # map!(x->Int(x/fact), values(res))
    # return res
    tmpB = zeros(Int, m)
    tmpC = zeros(Int, m)
    cycle = zeros(Int, m)
    # cycles = [zeros(Int, m) for i = 1:length(res)]
    # tmpM = BitMatrix(zeros(Bool, m, m))
    # M = BitArray(undef, m, m)
    


    labelledRes = Dict{Vector{Int8},Int}()    
    @showprogress 1 "labeling canonically..." for (B, c) in res
        @views @inbounds cycle[B[1,:]] = B[2,:]
        labelCanonical!(cycle, tmpB, tmpC)
        cL = Int8.(cycle)
        labelledRes[cL] = get(labelledRes, cL, 0) + c
    end

    # unLabelledRes = Dict{Vector{Int8},Int}()
    # for (B, c) in res 
    #     @views @inbounds cycle[B[1,:]] = B[2,:]
    #     cL = Int8.(cycle)
    #     unLabelledRes[cL] = c
    # end

    
    # labelledRes = Dict{Vector{Int8},Int}()
    # while length(unLabelledRes) > 0
    #     (cycle, c) = pop!(unLabelledRes)
    #     or = orbit(cycle)
    #     labC = minimum(or)

    #     for cyc in or 
    #         if haskey(unLabelledRes, cyc)
    #             c += unLabelledRes[cyc]
    #             delete!(unLabelledRes, cyc)
    #         end
    #     end
    #     if c != 0
    #         labelledRes[labC] = c
    #     end
    # end
    
    

    filter!(x -> x.second != 0, labelledRes)

    # @show typeof(labelledRes)

    
    # display(labelledRes)

    return labelledRes
end

# ## incomplete idea
# function calcProductFullSym(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})#, T::type = Int64)

#     #TODO: remove fact (since it is always 1)

#     m = max(maximum(t1), maximum(t2))

#     lambda = t1.part
#     lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

#     tL = length(t1.part)

#     function pMat(p::Perm{T})
#         t = zeros(Int8, tL, tL)
#         for (i, d) in enumerate(p.d)
#             t[i, d] = 1
#         end
#         return t
#     end

#     function oDet(n::T)
#         return Dict(pMat(P) => sign(P) for P in SymmetricGroup(n))
#     end

#     function mink(As::Dict{Matrix{Int8},T}, Bs::Dict{Matrix{Int8},T})
#         t = Dict{Matrix{Int8},T}()
#         for (A, ca) in As
#             for (B, cb) in Bs
#                 C = A + B
#                 t[C] = get!(t, C, 0) + ca * cb
#             end
#         end
#         return t
#     end

#     PLambda = Dict{Matrix{Int8},T}(zeros(Int8, tL, tL) => 1)
#     for k = 1:m
#         t = lambda[k] - lambda[k+1]
#         if t > 0
#             oD = oDet(T(k))
#             map!(x -> x * factorial(k), values(oD))
#             for i = 1:t
#                 PLambda = mink(PLambda, oD)
#             end
#         end
#     end


#     limit = false
#     r(s, j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x -> x == s, t1[j, :])) : 0)
#     u(s, j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x -> x == s, t2[j, :])) : 0)


#     # factor out the size of the column stabilizer
#     fact = prod([factorial(t) for t in conj(t1.part).part])
#     map!(x -> T(x / fact), values(PLambda))

#     toVec(x) = vcat(vcat(x[1]...), vcat(x[2]...), vec(x[3]), vec(x[4]))


#     res = Dict()
#     for (B, c) in PLambda
#         res[([t1[i, t1[i, :].!=0] for i = 1:tL], [t2[i, t2[i, :].!=0] for i = 1:tL], B, zeros(Int8, m, m))] = c
#     end

#     return res

#     return PLambda


#     return labelledRes
# end

##
# ## SDP calc speed test
# using AbstractAlgebra
# k = 8
# lambda = rand(FlagSOS.biggerShapes(AbstractAlgebra.Partition(repeat([1],k))))

# @time K = FlagSOS.Kostka(lambda,AbstractAlgebra.Partition(repeat([1],k)))[2]

# t1 = rand(K)
# t2 = rand(K)

# # @profview res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false)

# @time for i = 1:10; res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false); end
# # t = @timed res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false)

# # @show t.time

# ## 30-33s
# k = 4
# @time for i = 1:100
#     lambda = rand(FlagSOS.biggerShapes(AbstractAlgebra.Partition(repeat([1],k))))
#     @show i
#     K = FlagSOS.Kostka(lambda,AbstractAlgebra.Partition(repeat([1],k)))[2]

#     t1 = rand(K)
#     t2 = rand(K)
#     res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false)
# end

# ## 39s
# k = 7
# @time for i = 1:1000
#     lambda = rand(FlagSOS.biggerShapes(AbstractAlgebra.Partition(repeat([1],k))))
#     @show i
#      K = FlagSOS.Kostka(lambda,AbstractAlgebra.Partition(repeat([1],k)))[2]

#     t1 = rand(K)
#     t2 = rand(K)
#     res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false)
# end


# ## Fast sym polytab prod tests 
# using Combinatorics, AbstractAlgebra, SparseArrays, LinearAlgebra


# ##
# k = 9

# mu = AbstractAlgebra.Partition(repeat([1],k))
# mu = AbstractAlgebra.Partition([4,3,3,2,1,1])

# # mu = AbstractAlgebra.Partition([6,4,3,2,1,1,1])

# # for t = 1:10

# lambda = rand(FlagSOS.biggerShapes(mu))

# K = FlagSOS.Kostka(lambda,mu)[2]

# t1 = rand(K)
# t2 = rand(K)

# display(t1)
# display(t2)

# # @time for i = 1:10; res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false); end

# # @time for i = 1:10; res2 = test2(t1, t2); end

# # @time res = FlagSOS.symPolytabloidProduct(t1,t2,t1.part.part,false)

# # @time res2 = test2(t1, t2)
# # GC.gc()

# # @time res3 = test3(t1, t2)


# GC.gc()
# # both dir combined
# @time res5 = test5(t1, t2);

# GC.gc()
# # first one dir, then the other
# @time res6 = test6(t1, t2);

# GC.gc()
# # first one dir, then the other
# # s loop reversed
# # fastest? (on average)
# @time test7(t1, t2);

# # with partial symetry exploitation
# GC.gc()
# # @time test8(t1, t2)

# # @show Base.summarysize.([res2, res3, res4])

# # @show Base.summarysize(res3)
# ##
# pLi = 5

# @time res8 = test8(t1, t2)
# # @time res8b = test8b(t1, t2)
# # @time res8c = test8c(t1, t2)

# # @show length.([res8, res8b, res8c])

# #
# # symmetrized res8
# res8Sym = Dict{Matrix{Int8}, Int64}()
# for (ov, c) in res8
#     ovSym = reshape(maximum(vec(circshift(ov,(i,j))) for i = 1:k, j = 1:k), size(ov))
#     res8Sym[ovSym] = get!(res8Sym, ovSym, 0) + c
# end
# filter!(x->x[2]!=0, res8Sym)

# # display(res8Sym)
# @show length(res8Sym)

# #
# GC.gc()
# # @time res9 = test9(t2,t1)
# GC.gc()
# @time res9b = test9b(t2,t1)

# # display(res9)

# # if length(res8Sym) != length(res9)
# #     println("WRONG SIZE")
# # else
# #     @show unique(res8Sym[k]/res9[k] for k in keys(res8Sym))
# # end


# @show res8Sym == res9b
# # if res8Sym != res9
# #     break 
# # end
# # end
# ##




# function test8b(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})

#     m = max(maximum(t1), maximum(t2))


#     lambda = t1.part
#     lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

#     function move(A::Matrix{Int8},i::Int,j::Int,k::Int,l::Int)
#         res = copy(A)
#         res[i,j] -= 1
#         res[k,l] += 1
#         return res
#     end

#     function pMat(p::Perm{Int64})
#         t = zeros(Int8,m,m)
#         for (i,d) in enumerate(p.d)
#             t[i,d] = 1
#         end 
#         return t
#     end

#     function oDet(n::Int)
#         return Dict(pMat(P)=>sign(P) for P in SymmetricGroup(n))
#     end

#     function mink(As::Dict{Matrix{Int8}, Int64},Bs::Dict{Matrix{Int8}, Int64})
#         t = Dict{Matrix{Int8}, Int64}()
#         for (A,ca) in As
#             for (B, cb) in Bs
#                 C = A+B
#                 t[C] = get!(t, C, 0) + ca*cb
#             end
#         end
#         return t
#     end

#     res = Dict{Matrix{Int8}, Int64}(zeros(Int8,m,m) => 1)
#     for k = 1:m
#         t = lambda[k] - lambda[k+1]
#         if t > 0
#             oD = oDet(k)
#             map!(x->x*factorial(k), values(oD))
#             for i = 1:t
#                 res = mink(res, oD)
#             end
#         end
#     end

#     # res = Dict(collect(res)[pLi])
#     # display(res)

#     limit = false
#     r(s,j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x->x == s, t1[j,:])) : 0)
#     u(s,j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x->x == s, t2[j,:])) : 0)


#     fact = 1
#     for j = m-1:-1:1
#         # @show j
#         # First loop depends on order, second does not matter
#         for s = m:-1:j+1
#             # @show s
#             # @show (s,j)

#             usj = u(s,j)

#             fact *= factorial(usj)

#             for k = 1:usj
#                 res2 = Dict{Matrix{Int8}, Int64}()
#                 for (B,c) in res
#                     for nz in findall(!iszero, B[:,j])
#                         B2 = move(B, nz,j,nz,s)
#                         res2[B2] = get!(res2, B2, 0) + c*B[nz, j]
#                     end
#                 end
#                 res = res2
#             end
#         end
#     end


#     res2 = Dict{Matrix{Int8}, Int64}()
#     for (B,c) in res
#         B2 = reshape(minimum([vec(circshift(B,(0,i))) for i = 0:size(B,1)-1]), size(B))
#         res2[B2] = get!(res2, B2, 0) + c
#     end
#     filter!(x->x[2] != 0, res2)
#     res = res2


#     for j = m-1:-1:1
#         # @show j
#         # First loop depends on order, second does not matter
#         for s = m:-1:j+1
#             # @show s
#             # @show (s,j)

#             rsj = r(s,j)

#             fact *= factorial(rsj)


#             # res = unique(res)

#             for k = 1:rsj
#                 res2 = Dict{Matrix{Int8}, Int64}()
#                 for (B,c) in res
#                     for nz in findall(!iszero, B[j,:])
#                         B2 = move(B, j,nz,s,nz)
#                         res2[B2] = get!(res2, B2, 0) + c*B[j,nz]
#                     end
#                 end
#                 res = res2
#             end
#             res2 = Dict{Matrix{Int8}, Int64}()
#             for (B,c) in res
#                 B2 = reshape(minimum([vec(circshift(B,(0,i))) for i = 0:size(B,1)-1]), size(B))
#                 res2[B2] = get!(res2, B2, 0) + c
#             end
#             filter!(x->x[2] != 0, res2)
#             res = res2
#         end
#     end

#     map!(x->Int(x/fact), values(res))

#     return res
# end

# function test8c(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})

#     m = max(maximum(t1), maximum(t2))


#     lambda = t1.part
#     lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

#     function move(A::Matrix{Int8},i::Int,j::Int,k::Int,l::Int)
#         res = copy(A)
#         res[i,j] -= 1
#         res[k,l] += 1
#         return res
#     end

#     function pMat(p::Perm{Int64})
#         t = zeros(Int8,m,m)
#         for (i,d) in enumerate(p.d)
#             t[i,d] = 1
#         end 
#         return t
#     end

#     function oDet(n::Int)
#         return Dict(pMat(P)=>sign(P) for P in SymmetricGroup(n))
#     end

#     function mink(As::Dict{Matrix{Int8}, Int64},Bs::Dict{Matrix{Int8}, Int64})
#         t = Dict{Matrix{Int8}, Int64}()
#         for (A,ca) in As
#             for (B, cb) in Bs
#                 C = A+B
#                 t[C] = get!(t, C, 0) + ca*cb
#             end
#         end
#         return t
#     end

#     res = Dict{Matrix{Int8}, Int64}(zeros(Int8,m,m) => 1)
#     for k = 1:m
#         t = lambda[k] - lambda[k+1]
#         if t > 0
#             oD = oDet(k)
#             map!(x->x*factorial(k), values(oD))
#             for i = 1:t
#                 res = mink(res, oD)
#             end
#         end
#     end

#     # res = Dict(collect(res)[pLi])
#     # display(res)

#     limit = false
#     r(s,j) = ((j <= length(t1.part) && (!limit || j > 1)) ? (count(x->x == s, t1[j,:])) : 0)
#     u(s,j) = ((j <= length(t2.part) && (!limit || j > 1)) ? (count(x->x == s, t2[j,:])) : 0)


#     fact = 1
#     for j = m-1:-1:1
#         # @show j
#         # First loop depends on order, second does not matter
#         for s = m:-1:j+1
#             # @show s
#             # @show (s,j)

#             usj = u(s,j)

#             fact *= factorial(usj)

#             for k = 1:usj
#                 res2 = Dict{Matrix{Int8}, Int64}()
#                 for (B,c) in res
#                     for nz in findall(!iszero, B[:,j])
#                         B2 = move(B, nz,j,nz,s)
#                         res2[B2] = get!(res2, B2, 0) + c*B[nz, j]
#                     end
#                 end
#                 res = res2
#             end
#         end
#     end


#     res2 = Dict{Matrix{Int8}, Int64}()
#     for (B,c) in res
#         B2 = reshape(minimum([vec(circshift(B,(0,i))) for i = 0:size(B,1)-1]), size(B))
#         res2[B2] = get!(res2, B2, 0) + c
#     end
#     filter!(x->x[2] != 0, res2)
#     res = res2


#     for j = m-1:-1:1
#         # @show j
#         # First loop depends on order, second does not matter
#         for s = m:-1:j+1
#             # @show s
#             # @show (s,j)

#             rsj = r(s,j)

#             fact *= factorial(rsj)


#             # res = unique(res)

#             for k = 1:rsj
#                 res2 = Dict{Matrix{Int8}, Int64}()
#                 for (B,c) in res
#                     for nz in findall(!iszero, B[j,:])
#                         B2 = move(B, j,nz,s,nz)
#                         res2[B2] = get!(res2, B2, 0) + c*B[j,nz]
#                     end
#                 end
#                 res = res2
#             end

#         end
#         res2 = Dict{Matrix{Int8}, Int64}()
#         for (B,c) in res
#             B2 = reshape(minimum([vec(circshift(B,(0,i))) for i = 0:size(B,1)-1]), size(B))
#             res2[B2] = get!(res2, B2, 0) + c
#         end
#         filter!(x->x[2] != 0, res2)
#         res = res2
#     end

#     map!(x->Int(x/fact), values(res))

#     return res
# end



# # exploiting symmetry fully every step, only for Cm sym
# function test9(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})

#     m = max(maximum(t1), maximum(t2))


#     lambda = t1.part
#     lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

#     # function move(A::Matrix{Int8},i::Int,j::Int,k::Int,l::Int)
#     #     res = copy(A)
#     #     res[i,j] -= 1
#     #     res[k,l] += 1
#     #     return res
#     # end

#     function pMat(p::Perm{Int64})
#         t = zeros(Int8,m,m)
#         for (i,d) in enumerate(p.d)
#             t[i,d] = 1
#         end 
#         return t
#     end

#     function oDet(n::Int)
#         return Dict(pMat(P)=>sign(P) for P in SymmetricGroup(n))
#     end

#     function mink(As::Dict{Matrix{Int8}, Int64},Bs::Dict{Matrix{Int8}, Int64})
#         t = Dict{Matrix{Int8}, Int64}()
#         for (A,ca) in As
#             for (B, cb) in Bs
#                 C = A+B
#                 t[C] = get!(t, C, 0) + ca*cb
#             end
#         end
#         return t
#     end

#     res = Dict{Matrix{Int8}, Int64}(zeros(Int8,m,m) => 1)
#     for k = 1:m
#         t = lambda[k] - lambda[k+1]
#         if t > 0
#             oD = oDet(k)
#             map!(x->x*factorial(k), values(oD))
#             for i = 1:t
#                 res = mink(res, oD)
#             end
#         end
#     end


#     function translate(ov)
#         # display( ov)
#         # ov = ov'
#         m1 = zeros(Int8, m, m)
#         m1[diagind(m1)] .= 1
#         for i = 1:size(t1, 1)
#             if length(t1[i,:]) > 1
#                 for (a,b) in combinations(t1[i,:],2)
#                     if a!= 0 && b != 0
#                         m1[a,b] = 1
#                         m1[b,a] = 1
#                     end
#                 end
#             end
#         end

#         # @show(m1)


#         covered = [0 for i = 1:size(t1,1)]
#         permGroups = [[] for i = 1:size(t1,1)]
#         for i = 1:size(t2, 1)
#             for j = 1:size(t2, 1)
#                 for k = 1:ov[i,j]
#                     covered[i]+=1
#                     push!(permGroups[j], t1[i, covered[i]])
#                 end
#             end
#         end

#         # @show covered
#         # @show permGroups

#         m2 = zeros(Int8, m, m)
#         m2[diagind(m2)] .= 1
#         # for c in permGroups
#         #     if length(c) > 1
#         #         for (a,b) in combinations(c,2)
#         #             m2[a,b] = 1
#         #             m2[b,a] = 1
#         #         end
#         #     end
#         # end
#         for i = 1:size(t2, 1)
#             if length(t2[i,:]) > 1
#                 for (a,b) in combinations(t2[i,:],2)
#                     if a!= 0 && b != 0
#                         m2[a,b] = 1
#                         m2[b,a] = 1
#                     end
#                 end
#             end
#         end

#         # @show m2
#         # @show vec(t1')
#         # @show permGroups

#         m3 = zeros(Int8, m, m)
#         # @show conj(t1).fill
#         # @show vcat(permGroups...)

#         for (i,c) in enumerate(t2.fill)
#             m3[c,vcat(permGroups...)[i]] = 1
#         end

#         # @show permGroups
#         # @show m1
#         # @show m2
#         # @show m3

#         return (m1, m2, m3)
#     end

#     # countLabel = 0

#     function labelCan(tOv)
#         # @show countLabel += 1
#         # @show tOv
#         @views minV = maximum(vcat(vec(circshift(tOv[1], (i,i))),vec(circshift(tOv[2], (j,j))),vec(circshift(tOv[3], (j,i)))) for i = 1:m, j = 1:m)
#         return (reshape(minV[1:m^2],m,m), reshape(minV[m^2+1:2*m^2],m,m),reshape(minV[2*m^2+1:3*m^2],m,m))
#     end

#     # display(res)

#     res2 = Dict([])
#     for (i,(ov, c)) in enumerate(res)
#         # @show ov
#         # @show c
#         if true#i == pLi
#             # @show ov
#             # @show translate(ov)
#             tOv = labelCan(translate(ov))
#             res2[tOv] = get!(res2, tOv, 0) + c
#         end
#     end
#     res = res2

#     # display(res)

#     function labelRes(res)
#         res2 = Dict([])
#         for (ov, c) in res
#             tOv = labelCan(ov)
#             res2[tOv] = get!(res2, tOv, 0) + c
#         end
#         # @show length(res2)
#         @show (length(res), length(res2))
#         return res2
#     end

#     # for i = 1:sum(j-1 for j in t1.part)

#     #     # direction 1
#     #     res2 = Dict()

#     #     @show length(res)

#     #     for (ov, c) in res
#     #         # println("Split ov1: $ov, $c")
#     #         ovNoDiag = copy(ov[1])
#     #         ovNoDiag[diagind(ovNoDiag)] .= 0
#     #         pos = findfirst(ovNoDiag .== 1)[1]

#     #         allPos = findall(ov[1][pos,:] .== 1)
#     #         j=allPos[1]
#     #         newOv1 = copy(ov[1])
#     #         newOv2 = copy(ov[2])


#     #         otherPos = setdiff(allPos,[j])
#     #         newOv1[j,otherPos] .= 0
#     #         newOv1[otherPos,j] .= 0

#     #         for k in allPos
#     #             newOv3 = copy(ov[3])
#     #             newOv3[:,k] = ov[3][:,j]
#     #             newOv3[:,j] = ov[3][:,k]

#     #             # @show newOv1

#     #             # @show newOv2

#     #             # @show newOv3

#     #             res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#     #         end
#     #     end


#     #     res = res2
#     #     # res = labelRes(res2)
#     #     @show length(res)

#     #     # direction 2
#     #     res2 = Dict()
#     #     for (ov, c) in res
#     #         # println("Split ov2: $ov, $c")
#     #         ovNoDiag = copy(ov[2])
#     #         ovNoDiag[diagind(ovNoDiag)] .= 0
#     #         pos = findfirst(ovNoDiag .== 1)[1]

#     #         allPos = findall(ov[2][pos,:] .== 1)
#     #         j=allPos[1]
#     #         newOv1 = copy(ov[1])
#     #         newOv2 = copy(ov[2])


#     #         otherPos = setdiff(allPos,[j])
#     #         newOv2[j,otherPos] .= 0
#     #         newOv2[otherPos,j] .= 0

#     #         for k in allPos
#     #             newOv3 = copy(ov[3])
#     #             newOv3[k,:] = ov[3][j,:]
#     #             newOv3[j,:] = ov[3][k,:]


#     #             # @show newOv1

#     #             # @show newOv2

#     #             # @show newOv3

#     #             res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#     #         end
#     #     end
#     #     res = labelRes(res2)
#     #     @show length(res)

#     # end

#     for i = 1:sum(j-1 for j in t1.part)

#         # direction 1
#         res2 = Dict()

#         # @show length(res)

#         for (ov, c) in res
#             # println("Split ov1: $ov, $c")
#             ovNoDiag = copy(ov[1])
#             ovNoDiag[diagind(ovNoDiag)] .= 0
#             pos = findfirst(ovNoDiag .== 1)[1]

#             allPos = findall(ov[1][pos,:] .== 1)
#             j=allPos[1]
#             newOv1 = copy(ov[1])
#             newOv2 = copy(ov[2])


#             otherPos = setdiff(allPos,[j])
#             newOv1[j,otherPos] .= 0
#             newOv1[otherPos,j] .= 0

#             for k in allPos
#                 newOv3 = copy(ov[3])
#                 newOv3[:,k] = ov[3][:,j]
#                 newOv3[:,j] = ov[3][:,k]

#                 # @show newOv1

#                 # @show newOv2

#                 # @show newOv3

#                 res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#             end
#         end

#         res = labelRes(res2)
#         # @show length(res)

#     end

#     for i = 1:sum(j-1 for j in t1.part)



#         # direction 2
#         res2 = Dict()
#         for (ov, c) in res
#             # println("Split ov2: $ov, $c")
#             ovNoDiag = copy(ov[2])
#             ovNoDiag[diagind(ovNoDiag)] .= 0
#             pos = findfirst(ovNoDiag .== 1)[1]

#             allPos = findall(ov[2][pos,:] .== 1)
#             j=allPos[1]
#             newOv1 = copy(ov[1])
#             newOv2 = copy(ov[2])


#             otherPos = setdiff(allPos,[j])
#             newOv2[j,otherPos] .= 0
#             newOv2[otherPos,j] .= 0

#             for k in allPos
#                 newOv3 = copy(ov[3])
#                 newOv3[k,:] = ov[3][j,:]
#                 newOv3[j,:] = ov[3][k,:]


#                 # @show newOv1

#                 # @show newOv2

#                 # @show newOv3

#                 res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#             end
#         end
#         res = labelRes(res2)
#         # @show length(res)

#     end



#     return Dict(ov[3]=>c for (ov, c) in res if c != 0)
# end

# # even more symmetry
# function test9b(t1::AbstractAlgebra.Generic.YoungTableau{Int64}, t2::AbstractAlgebra.Generic.YoungTableau{Int64})

#     m = max(maximum(t1), maximum(t2))


#     lambda = t1.part
#     lambda = vcat(lambda, [zero(lambda[1]) for i = length(lambda)+1:m+1])

#     # function move(A::Matrix{Int8},i::Int,j::Int,k::Int,l::Int)
#     #     res = copy(A)
#     #     res[i,j] -= 1
#     #     res[k,l] += 1
#     #     return res
#     # end

#     function pMat(p::Perm{Int64})
#         t = zeros(Int8,m,m)
#         for (i,d) in enumerate(p.d)
#             t[i,d] = 1
#         end 
#         return t
#     end

#     function oDet(n::Int)
#         return Dict(pMat(P)=>sign(P) for P in SymmetricGroup(n))
#     end

#     function mink(As::Dict{Matrix{Int8}, Int64},Bs::Dict{Matrix{Int8}, Int64})
#         t = Dict{Matrix{Int8}, Int64}()
#         for (A,ca) in As
#             for (B, cb) in Bs
#                 C = A+B
#                 t[C] = get!(t, C, 0) + ca*cb
#             end
#         end
#         return t
#     end

#     res = Dict{Matrix{Int8}, Int64}(zeros(Int8,m,m) => 1)
#     for k = 1:m
#         t = lambda[k] - lambda[k+1]
#         if t > 0
#             oD = oDet(k)
#             map!(x->x*factorial(k), values(oD))
#             for i = 1:t
#                 res = mink(res, oD)
#             end
#         end
#     end


#     function translate(ov)
#         # display( ov)
#         # ov = ov'
#         m1 = zeros(Int8, m, m)
#         m1[diagind(m1)] .= 1
#         for i = 1:size(t1, 1)
#             if length(t1[i,:]) > 1
#                 for (a,b) in combinations(t1[i,:],2)
#                     if a!= 0 && b != 0
#                         m1[a,b] = 1
#                         m1[b,a] = 1
#                     end
#                 end
#             end
#         end

#         # @show(m1)


#         covered = [0 for i = 1:size(t1,1)]
#         permGroups = [[] for i = 1:size(t1,1)]
#         for i = 1:size(t2, 1)
#             for j = 1:size(t2, 1)
#                 for k = 1:ov[i,j]
#                     covered[i]+=1
#                     push!(permGroups[j], t1[i, covered[i]])
#                 end
#             end
#         end

#         # @show covered
#         # @show permGroups

#         m2 = zeros(Int8, m, m)
#         m2[diagind(m2)] .= 1
#         # for c in permGroups
#         #     if length(c) > 1
#         #         for (a,b) in combinations(c,2)
#         #             m2[a,b] = 1
#         #             m2[b,a] = 1
#         #         end
#         #     end
#         # end
#         for i = 1:size(t2, 1)
#             if length(t2[i,:]) > 1
#                 for (a,b) in combinations(t2[i,:],2)
#                     if a!= 0 && b != 0
#                         m2[a,b] = 1
#                         m2[b,a] = 1
#                     end
#                 end
#             end
#         end

#         # @show m2
#         # @show vec(t1')
#         # @show permGroups

#         m3 = zeros(Int8, m, m)
#         # @show conj(t1).fill
#         # @show vcat(permGroups...)

#         for (i,c) in enumerate(t2.fill)
#             m3[c,vcat(permGroups...)[i]] = 1
#         end

#         # @show permGroups
#         # @show m1
#         # @show m2
#         # @show m3

#         return (m1, m2, m3)
#     end

#     # countLabel = 0

#     function labelCan(tOv)

#         maxM1 = maximum(vec(circshift(tOv[1], (i,i))) for i = 1:m)
#         maxM1Ind = findall((x->x==maxM1).([vec(circshift(tOv[1], (i,i))) for i = 1:m]))
#         maxM1Mat = reshape(maxM1, m,m)
#         grps1 = []
#         cur = 1
#         while cur <= m
#             if length(grps1) == 0 || all(!in(cur, g) for g in grps1)
#                 push!(grps1, findall(maxM1Mat[cur,:] .== 1))
#             end
#             cur += 1
#             # cur += length(grps1[end])
#         end

#         if length(unique(vcat(grps1...))) != m
#             @error("$grps1")
#         end

#         # @show maxM1Mat
#         # @show grps1
#         sym1 = []
#         for p in Iterators.product([[maxM1Ind]; [SymmetricGroup(length(I)) for I in grps1]]...)
#             # @show p
#             per = perm(1:m)
#             for (i, I) in enumerate(grps1)
#                 ptmp = collect(1:m)
#                 for (j, c) in enumerate(p[i+1].d)
#                     ptmp[I[j]] = I[c]
#                 end
#                 per *= perm(ptmp)
#             end
#             push!(sym1, per*(perm([2:m; 1])^(m-p[1])))

#         end
#         # display(sym1)


#         maxM2 = maximum(vec(circshift(tOv[2], (i,i))) for i = 1:m)
#         maxM2Ind = findall((x->x==maxM2).([vec(circshift(tOv[2], (i,i))) for i = 1:m]))
#         maxM2Mat = reshape(maxM2, m,m)
#         grps2 = []
#         cur = 1
#         while cur <= m
#             if length(grps2) == 0 || all(!in(cur, g) for g in grps2)
#                 push!(grps2, findall(maxM2Mat[cur,:] .== 1))
#             end
#             cur += 1
#             # cur += length(grps1[end])
#         end

#         if length(unique(vcat(grps2...))) != m
#             @error("$grps2")
#         end

#         sym2 = []
#         for p in Iterators.product([[maxM2Ind]; [SymmetricGroup(length(I)) for I in grps2]]...)
#             # @show p
#             per = perm(1:m)
#             for (i, I) in enumerate(grps2)
#                 ptmp = collect(1:m)
#                 for (j, c) in enumerate(p[i+1].d)
#                     ptmp[I[j]] = I[c]
#                 end
#                 per *= perm(ptmp)
#             end
#             push!(sym2, per*(perm([2:m; 1])^(m-p[1])))

#         end

#         maxM3 = reshape(maximum(vec(tOv[3][q.d,p.d]) for p in sym1, q in sym2),m,m)

#         # display(maxM1Mat)

#         return (maxM1Mat, maxM2Mat, maxM3)
#     end

#     # display(res)

#     res2 = Dict([])
#     for (i,(ov, c)) in enumerate(res)
#         # @show ov
#         # @show c
#         if true#i == pLi
#             # @show ov
#             # @show translate(ov)
#             tOv = labelCan(translate(ov))
#             res2[tOv] = get!(res2, tOv, 0) + c
#         end
#     end
#     res = res2

#     display(length(res))

#     function labelRes(res)
#         res2 = Dict([])
#         for (ov, c) in res
#             tOv = labelCan(ov)
#             res2[tOv] = get!(res2, tOv, 0) + c
#         end
#         # @show length(res2)
#         @show (length(res), length(res2))
#         # display(rand(res2)[1][1])
#         # display(rand(res2)[1][2])

#         # rO = rand(res2)[1]
#         # display.(rO)
#         return res2
#     end

#     for i = 1:sum(j-1 for j in t1.part)

#         # direction 1
#         res2 = Dict()

#         # @show length(res)

#         for (ov, c) in res
#             # println("Split ov1: $ov, $c")
#             ovNoDiag = copy(ov[1])
#             ovNoDiag[diagind(ovNoDiag)] .= 0
#             pos = findfirst(ovNoDiag .== 1)[1]

#             allPos = findall(ov[1][pos,:] .== 1)
#             j=allPos[1]
#             newOv1 = copy(ov[1])
#             newOv2 = copy(ov[2])


#             otherPos = setdiff(allPos,[j])
#             newOv1[j,otherPos] .= 0
#             newOv1[otherPos,j] .= 0

#             for k in allPos
#                 newOv3 = copy(ov[3])
#                 newOv3[:,k] = ov[3][:,j]
#                 newOv3[:,j] = ov[3][:,k]

#                 # @show newOv1

#                 # @show newOv2

#                 # @show newOv3

#                 res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#             end
#         end


#         # res = res2
#         res = labelRes(res2)
#         # @show length(res)

#         # # direction 2
#         # res2 = Dict()
#         # for (ov, c) in res
#         #     # println("Split ov2: $ov, $c")
#         #     ovNoDiag = copy(ov[2])
#         #     ovNoDiag[diagind(ovNoDiag)] .= 0
#         #     pos = findfirst(ovNoDiag .== 1)[1]

#         #     allPos = findall(ov[2][pos,:] .== 1)
#         #     j=allPos[1]
#         #     newOv1 = copy(ov[1])
#         #     newOv2 = copy(ov[2])


#         #     otherPos = setdiff(allPos,[j])
#         #     newOv2[j,otherPos] .= 0
#         #     newOv2[otherPos,j] .= 0

#         #     for k in allPos
#         #         newOv3 = copy(ov[3])
#         #         newOv3[k,:] = ov[3][j,:]
#         #         newOv3[j,:] = ov[3][k,:]


#         #         # @show newOv1

#         #         # @show newOv2

#         #         # @show newOv3

#         #         res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#         #     end
#         # end
#         # res = labelRes(res2)
#         # # @show length(res)

#     end

#     for i = 1:sum(j-1 for j in t1.part)

#         # # direction 1
#         # res2 = Dict()

#         # # @show length(res)

#         # for (ov, c) in res
#         #     # println("Split ov1: $ov, $c")
#         #     ovNoDiag = copy(ov[1])
#         #     ovNoDiag[diagind(ovNoDiag)] .= 0
#         #     pos = findfirst(ovNoDiag .== 1)[1]

#         #     allPos = findall(ov[1][pos,:] .== 1)
#         #     j=allPos[1]
#         #     newOv1 = copy(ov[1])
#         #     newOv2 = copy(ov[2])


#         #     otherPos = setdiff(allPos,[j])
#         #     newOv1[j,otherPos] .= 0
#         #     newOv1[otherPos,j] .= 0

#         #     for k in allPos
#         #         newOv3 = copy(ov[3])
#         #         newOv3[:,k] = ov[3][:,j]
#         #         newOv3[:,j] = ov[3][:,k]

#         #         # @show newOv1

#         #         # @show newOv2

#         #         # @show newOv3

#         #         res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#         #     end
#         # end


#         # # res = res2
#         # res = labelRes(res2)
#         # @show length(res)

#         # direction 2
#         res2 = Dict()
#         for (ov, c) in res
#             # println("Split ov2: $ov, $c")
#             ovNoDiag = copy(ov[2])
#             ovNoDiag[diagind(ovNoDiag)] .= 0
#             pos = findfirst(ovNoDiag .== 1)[1]

#             allPos = findall(ov[2][pos,:] .== 1)
#             j=allPos[1]
#             newOv1 = copy(ov[1])
#             newOv2 = copy(ov[2])


#             otherPos = setdiff(allPos,[j])
#             newOv2[j,otherPos] .= 0
#             newOv2[otherPos,j] .= 0

#             for k in allPos
#                 newOv3 = copy(ov[3])
#                 newOv3[k,:] = ov[3][j,:]
#                 newOv3[j,:] = ov[3][k,:]


#                 # @show newOv1

#                 # @show newOv2

#                 # @show newOv3

#                 res2[(newOv1, newOv2, newOv3)] = get!(res2, (newOv1, newOv2, newOv3), 0) + c
#             end
#         end
#         res = labelRes(res2)
#         # @show length(res)

#     end



#     return Dict(ov[3]=>c for (ov, c) in res if c != 0)
# end


# ## sparse matrix approach: Try to write all steps of the algorithm as sparse matrix - sparse vector multiplications. Not finished

# # all overlaps that may appear in the result. Does miss intermediate matrices.
# # function allOverlaps(t1, t2)
# #     m = max(maximum(t1), maximum(t2))
# #     fill1 = [count(t1.==i) for i = 1:m]
# #     fill2 = [count(t2.==i) for i = 1:m]

# #     lastStep = [zeros(Int, m, m)]


# #     for i = 1:m
# #         nextStep = []
# #         for c in combinations(1:m+fill1[i]-1, m-1)
# #             d = [0,c...,m+fill1[i]]
# #             newRow = [d[i+1] - d[i] - 1 for i = 1:m]
# #             for A in lastStep
# #                 newA = copy(A)
# #                 newA[i,:] = newRow
# #                 if all(vec(sum(newA, dims=1)) .<= fill2)
# #                     push!(nextStep, newA)
# #                 end
# #             end
# #         end
# #         lastStep = nextStep
# #     end
# #     return lastStep
# # end

# # ov = allOverlaps(t1,t2)

# # filter!(x->x[1,1] == 1, ov)





# # unique(reshape(maximum([vec(circshift(o,(i,j))) for i = 1:k, j=1:k]),k,k) for o in ov)

## Comparison to Sven's results
# m = 9
# basis = decomposeModule(m)
# # P = Partition([6,3])
# # svenRes = CrossingBlockDiagOnly(m)

# using FileIO
# # save("data/svenBlockDiag9.jld2", "blks", svenRes)
# svenRes = load("data/svenBlockDiag9.jld2")["blks"]

# myRes = blockDiagonalize(basis)

# # blockSize = length(basis[(P,-1)])

# ##
# # blockPPlus = svenRes[2][1]
# svenTranslate = Dict(c=>labelCanonical(k) for (k,c) in svenRes[1])


# svenDataTranslated = Dict()

# for (i, P) in enumerate(keys(basis))
#     curBlock = svenRes[2][i]
#     svenDataTranslated[P] = Dict()
#     blockSize = length(basis[P])

#     for (num, c) in svenTranslate
#         svenDataTranslated[P][c] = zeros(Int,blockSize, blockSize)
#         for i = 1:blockSize
#             for j = i:blockSize
#                 if haskey(curBlock[[i,j]], num)
#                     svenDataTranslated[P][c][i,j] = Int(0.5*curBlock[[i,j]][num])
#                 end
#             end
#         end
#     end
#     filter!(p->!iszero(p.second), svenDataTranslated[P])

# end


# # A = myRes[(P,1)]

# test = [unique(vcat([unique(myRes[P][c] ./svenDataTranslated[P][c])   for c in keys(myRes[P])]...)) for P in keys(basis)]

# # c = rand(keys(A))

# # display(A[c])
# # display(svenDataTranslated[c])
# # display(A[c] ./ svenDataTranslated[c])

# myRes == svenDataTranslated



# ##
# for P in keys(myRes)
#     @show (length(keys(myRes[P])), length(keys(svenDataTranslated[P])))
#     if length(keys(myRes[P])) <  length(keys(svenDataTranslated[P]))
#         global testP = P 
#     end
#     for c in keys(myRes[P])
#         if myRes[P][c] != svenDataTranslated[P][c]
#             @show (P,c)
#         end
#     end
# end