## Code from FlagSOS.jl (TODO: use proper package once its public)

using AbstractAlgebra: YoungTableau, Partition, @perm_str, matrix_repr, AllParts
using AbstractAlgebra.Generic: collength, rowlength
using Combinatorics: combinations

function setTabEntry(Y,p,c)
    linPos = 0
    for i = 1:(p[1]-1)
        linPos += Y.part[i]
    end
    linPos += p[2]
    Y.fill[linPos] = c
end



# Calculating Kostka numbers
# number of semistandard tableaux of shape λ and content μ
# i.e. increasing rows, strictly increasing columns
# Also returns semistandard tableaux
function Kostka(λ,μ)
    if λ.n != μ.n
        return nothing
    end
    # Fill T up with num times c
    function fill(T,c, num)
        # Find potential rows with spaces
        potential = [0 for i=0:collength(T, 1,1)]
        lastNZ = [0 for i=0:collength(T, 1,1)]
        for i = 1:length(potential)
            # Find first zero in row i of T
            if T[i, 1] != 0
                lastNZ[i] = 1

                while rowlength(T, i,lastNZ[i]) > 0 && T[i,lastNZ[i]] != 0
                    lastNZ[i] += 1
                end
                if T[i, lastNZ[i]] == 0
                    lastNZ[i] -= 1
                end
                if  rowlength(T, i,lastNZ[i]) == 0
                    continue
                end
            end

            potential[i] = min(rowlength(T, i,1) + 1, (i>1 ? lastNZ[i-1] : 999999999999)) - lastNZ[i]
        end
        # display((c, potential,lastNZ))

        if num > sum(potential .* (potential .> 0))
            return []
        end

        numPotentialRows = sum(potential .> 0)
        res = []

        function fillRows(TC, r, remNum)
            for i = max(0,remNum - sum(potential[r+1:end])):min(potential[r], remNum)
                Tres = YoungTableau(TC.part, copy(TC.fill))
                for col = (lastNZ[r]+1):(lastNZ[r]+i)
                    setTabEntry(Tres,(r, col), c)
                end
                if i == remNum
                    # display(Tres)
                    push!(res, Tres)
                else
                    fillRows(Tres, r+1, remNum - i)
                end
            end
        end

        fillRows(T, 1, num)
        return res
    end

    TLast = [YoungTableau(λ,[0 for i=1:λ.n])]
    for i = 1:length(μ.part)
        TCur = []
        for T in TLast
            TCur = union!(TCur, fill(T,i, μ[i]))
        end
        TLast = TCur
    end
    return (length(TLast), TLast)
end

#

function biggerShapes(λ)
    ## Find bigger shapes

    if length(λ) == 0
        return [λ]
    end

    firstRow = λ[1]
    res = [Partition([λ.n])]
    for i = 1:(λ.n - firstRow)
        for p in AllParts(i)
            if length(p) < length(λ) && p[1] <= λ.n - i
                tmp = vcat([λ.n-i], p)
                # display((tmp,i,p))
                valid = true
                if i == λ.n - firstRow
                    for j = 1:length(tmp)
                        if tmp[j] > λ[j]
                            break
                        end
                        if tmp[j] < λ[j]
                            valid = false
                            # display((tmp, λ))
                            break
                        end
                    end
                end
                if valid
                    push!(res, Partition(tmp))
                end
            end
        end
    end
    return sort!(res)
end

function semiTabPermMatV3(K, p)
    tmpRes = Dict()

    function applyPermToSemiTabV3(semiTab, p)

        t0 = YoungTableau(semiTab.part, [p[i] for i in semiTab.fill])
        # @show t0
    
        function isSemiUpToi(t, i)
            last = zeros(Int64, size(t,2))
    
            for j = 1:i
                current = vec(sum(t .== j, dims = 1))
                newLast = last + current
                if maximum(current) > 1 || !issorted(newLast, rev = true) 
                    return false
                end
                last = newLast
            end
            return true
        end
    
    
        function decomp(t, ignoreT)
            return get!(tmpRes, (t, ignoreT)) do
    
    
                iFixed = [t]
    
                for i = 1:maximum(t.fill)
                    tmp = []
                    for ti in iFixed
                        tiTmp = []
                        for j = 1:size(ti,1)
                            @inbounds @views posPositions = findall(ti[j,:] .>= i)
                            @inbounds @views numPos = count(ti[j,:] .== i)
                            rowOptions = []
                            for c in combinations(posPositions, numPos)
                                @inbounds @views tmpRow = copy(ti[j,:])
                                @inbounds tmpRow[c] .= i
                                @inbounds @views tmpRow[setdiff(posPositions, c)] = filter(x->x!=i, ti[j,posPositions])
                                push!(rowOptions, reshape(tmpRow, 1, :))
                            end
    
                            if j == 1
                                tiTmp = rowOptions
                            else
                                tiTmp = [vcat(a,b) for a in tiTmp for b in rowOptions]
                            end
                        end
                        tmp = vcat(tmp, filter(x->isSemiUpToi(x, i),tiTmp))
                    end
                    iFixed = tmp
                end
    
                # sort collumns and calculate signs
                res = Dict()
                for ti in iFixed
                    s = 1
                    for j = 1:size(ti,2)
                        col = ti[1:collength(t,1,j)+1,j]
                        if length(col) > 1
                            s *= (-1)^sum(c2 < c1 for (c1,c2) in combinations(col,2))
                        end
                        sCol = sort(col)
                        for i = 1:length(sCol)
                            ti[i,j] = sCol[i]
                        end
                    end
    
                    yti = YoungTableau(t.part, filter(x->x!=0,vec(ti')))
    
                    if !ignoreT || yti != t 
                        res[yti] = get(res, yti,0) + s
                    end
                end
                
                fullres = deepcopy(res)
    
                for (ti,c) in res
                    tmp = decomp(ti, true)
                    for (tj, c2) in tmp
                        fullres[tj] = get(fullres, tj, 0) - c*c2
                    end
                end
    
    
                fullres
            end
        end
        return collect(decomp(t0, false))
    end

    mat = zeros(Int64, length(K[2]), length(K[2]))
    println("")
    for (i, t) in enumerate(K[2])
        print("$i / $(K[1])\r")
        res = applyPermToSemiTabV3(t, p)
        if res !== nothing
            for r in res
                c = findfirst(x -> x == r[1], K[2])
                mat[c, i] += r[2]
            end
        end
    end
    return mat
end


function findColSpanningSetV3(A)

    F = qr(A, Val(true))

    r = sum(abs.(diag(F.R)) .> 0.00000001)

    if r == 0
        return []
    end

    return sort(F.p[1:r])

    # for sparse matrices?
    # return sort(qr(sparse(Float64.(A))).cpiv[1:rank(A)])

    

end


function generateGroup(generator, size, hasRepMat = false)
    if size == 1
        return generator
        # return [[]]
    end
    nextIndex = length(generator)
    if !hasRepMat
        res = copy(generator)
        lastIndex = 0

        while nextIndex < size
            union!(res, [f * g for f in res[lastIndex+1:length(res)], g in generator])
            lastIndex = nextIndex
            nextIndex = length(res)
            # display(nextIndex)
        end
        return res
    else
        representation = Dict([g[1] => g[2] for g in generator])
        res = [g[1] for g in generator]
        lastIndex = 0
        while nextIndex < size
            for f in res[lastIndex+1:length(res)]
                for g in generator
                    # p = g[1] * f
                    p = [f[g[1][i]] for i = 1:length(g[1])]
                    if !in(p, res)
                        push!(res, p)
                        merge!(
                            representation,
                            Dict([p => representation[f] * representation[g[1]]]),
                        )
                    end
                end
            end
            lastIndex = nextIndex
            nextIndex = length(res)
        end
        return (res, representation)
    end
end