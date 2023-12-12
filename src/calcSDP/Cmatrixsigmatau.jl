using LinearAlgebra
using Combinatorics
using JLD2, FileIO
using Printf
using DataStructures

setprecision(128)

function Csigmatau(m)
    Sigmavec = Array((1:m))
    Vertices = collect(permutations(Array((2:m))))
    for vertex in Vertices
        pushfirst!(vertex, 1)
    end
    DistanceDict = Dict()
    DistanceDict[Sigmavec] = 0
    Queue = PriorityQueue{Array{Int8},Int}()
    Queue[Sigmavec] = 0
    VisitedSet = Set{Array{Int8}}()

    while !isempty(Queue)
        #determine the current vertex
        CurrentDistance = peek(Queue)[2]
        currentvertex = dequeue!(Queue)
        push!(VisitedSet, currentvertex)

        #iterate over neighbors of current vertex, but only if they have not been visited.
        for i1 = 2:m-1
            NewNeighbor = deepcopy(currentvertex)
            NewNeighbor[[i1, i1 + 1]] = NewNeighbor[[i1 + 1, i1]]
            if !(NewNeighbor in VisitedSet)
                tempDistance = get(DistanceDict, NewNeighbor, CurrentDistance + 2)
                if tempDistance > CurrentDistance + 1
                    tempDistance = CurrentDistance + 1
                    DistanceDict[NewNeighbor] = tempDistance
                end
                Queue[NewNeighbor] = tempDistance
            end
        end
        #now swap 1 and 2
        NewNeighbor = vcat(Int8(1), currentvertex[3:end], currentvertex[2])
        if !(NewNeighbor in VisitedSet)
            tempDistance = get(DistanceDict, NewNeighbor, CurrentDistance + 2)
            if tempDistance > CurrentDistance + 1
                tempDistance = CurrentDistance + 1
                DistanceDict[NewNeighbor] = tempDistance
            end
            Queue[NewNeighbor] = tempDistance
        end
        #now swap m and 1
        NewNeighbor = vcat(Int8(1), currentvertex[m], currentvertex[2:end-1])
        if !(NewNeighbor in VisitedSet)
            tempDistance = get(DistanceDict, NewNeighbor, CurrentDistance + 2)
            if tempDistance > CurrentDistance + 1
                tempDistance = CurrentDistance + 1
                DistanceDict[NewNeighbor] = tempDistance
            end
            Queue[NewNeighbor] = tempDistance
        end
    end
    return DistanceDict
end


# label canonical without (A,B) -> (B,A) swap
labelDictDist = Dict{Vector{Int8},Vector{Int8}}()
# Optimize further! Currently faster with dict, probably possible to be faster without
function labelCanonicalDist(cycle::Vector{Int8}, best::Vector{Int8} = copy(cycle))
    # return get!(labelDictDist, cycle) do
    n = length(cycle)
    best .= cycle
    
    function lbl(c::Vector{Int8})
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
            for (a,b) in Iterators.zip(Iterators.flatten((i:n,1:i-1)), Iterators.flatten((1:(n-i+1),n-i+2:n)))
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
                @views @inbounds best[1:(n-i+1)] .= c[i:n]
                @views @inbounds best[n-i+2:end] .= c[1:i-1]
            end
        end
        return best
        # nothing
        # return minimum(((circshift(cycle,-i) .+ (n-cycle[i+1])) .% n) .+ 1 for i = 0:n-1)
    end

    # function lbl(cycle)
    #     return minimum(((circshift(cycle, -i) .+ (n - cycle[i+1])) .% n) .+ 1 for i = 0:n-1)
    # end

    # doubleInvCycle = n+1 .- reverse(cycle)

    lbl(cycle)
    reverse!(cycle)
    cycle .= n + 1 .- cycle
    lbl(cycle)
    cycle .= best
    return cycle

    # return min(lbl(cycle), lbl(doubleInvCycle))
    # end
end



function CsigmatauOrbits(m::Int)
    Sigmavec = Array{Int8}((1:m))
    # Vertices = collect(permutations(Array((2:m))));
    # for vertex in Vertices
    #     pushfirst!(vertex, 1);
    # end
    DistanceDict = Dict{Vector{Int8},Int}()
    DistanceDict[Sigmavec] = 0
    Queue = PriorityQueue{Array{Int8},Int}()
    Queue[Sigmavec] = 0
    VisitedSet = Set{Array{Int8}}()
    i = 1

    tmp = zeros(Int8, m)

    
    push!(VisitedSet, Sigmavec)

    while !isempty(Queue)
        #determine the current vertex
        CurrentDistance = peek(Queue)[2]
        currentvertex = dequeue!(Queue)
        # push!(VisitedSet, currentvertex)

        # println("Visited orbits: $(i)\r")
        i % 1000 == 1 && print("Visited orbits: $(i), Queue length: $(length(Queue))\r")
        i += 1
        # iterate over neighbors of current vertex, but only if they have not been visited.

        # More compact, but slower for some reason
        # swaps = [sort!([i, i == m ? 1 : i+1]) for i = 1:m]
        # indices = [vcat(1:a-1, b, a+1:b-1, a, b+1:m) for (a,b) in swaps]
        # NewNeighbors = Set{Vector{Int8}}([labelCanonical(currentvertex[s]) for s in indices])

        NewNeighbors = Set{Vector{Int8}}()

        for i1 = 2:m-1
            NewNeighbor = deepcopy(currentvertex)
            @inbounds @views NewNeighbor[[i1, i1 + 1]] .= NewNeighbor[[i1 + 1, i1]]
            push!(NewNeighbors, NewNeighbor)
        end
        #now swap 1 and 2
        @inbounds @views NewNeighbor = vcat(Int8(1), currentvertex[3:end], currentvertex[2])
        push!(NewNeighbors, NewNeighbor)

        #now swap m and 1
        @inbounds @views NewNeighbor = vcat(Int8(1), currentvertex[m], currentvertex[2:end-1])
        push!(NewNeighbors, NewNeighbor)


        # NewNeighborsInv = Set{Vector{Int8}}()
        # currentVerterxInv = vcat(Int8(1),currentvertex[end:-1:2])
        # for i1=2:m-1
        #     NewNeighbor = deepcopy(currentVerterxInv);
        #     NewNeighbor[[i1,i1+1]]=NewNeighbor[[i1+1,i1]]
        #     push!(NewNeighborsInv, NewNeighbor)
        # end
        # #now swap 1 and 2
        # NewNeighbor = vcat(Int8(1),currentVerterxInv[3:end],currentVerterxInv[2])
        # push!(NewNeighborsInv, NewNeighbor)

        # #now swap m and 1
        # NewNeighbor = vcat(Int8(1),currentVerterxInv[m],currentVerterxInv[2:end-1])
        # push!(NewNeighborsInv, NewNeighbor)

        NewNeighbors = unique(labelCanonicalDist(ne, tmp) for ne in NewNeighbors)
        # NewNeighborsInv = unique(labelCanonical.(NewNeighborsInv))

        # union!(NewNeighbors, NewNeighborsInv)
        # if NewNeighbors == NewNeighborsInv
        #     println("DIfferent neighbors!")
        #     display(currentvertex)
        #     display(NewNeighbors)
        #     display(NewNeighborsInv)
        # end

        setdiff!(NewNeighbors, VisitedSet)
        union!(VisitedSet, NewNeighbors)

        for NewNeighbor in NewNeighbors
            # tempDistance = get(DistanceDict, NewNeighbor, CurrentDistance + 2)
            # if tempDistance > CurrentDistance + 1
            #     tempDistance = CurrentDistance + 1
            #     DistanceDict[NewNeighbor] = tempDistance
            # end
            # Queue[NewNeighbor] = tempDistance
            DistanceDict[NewNeighbor] = CurrentDistance + 1
            Queue[NewNeighbor] = CurrentDistance + 1
        end
    end
    # @show length(VisitedSet)
    # return length(DistanceDict)

    M = BitArray(undef, m, m)
    tmp1 = zeros(Int, m)
    tmp2 = zeros(Int, m)
    tmp3 = zeros(Int, m)


    res = Dict{Vector{Int8}, Int}()
    for (k, v) in DistanceDict
        # M .= 0
        # for i = 1:m
        #     M[i, k[i]] = 1
        # end
        tmp1 = Int.(k)
        labelCanonical!(tmp1, tmp2, tmp3)
        res[Int8.(tmp1)] = v
    end
    return res
    # return Dict{Vector{Int8},Int}(labelCanonical(k) => v for (k, v) in DistanceDict)
    # return DistanceDict;
end


function CsigmatauOld(m)
    Sigmavec = Array((1:m))

    Vertices = collect(permutations(Array((2:m))))
    for vertex in Vertices
        pushfirst!(vertex, 1)
    end

    DistanceDict = Dict()
    DistanceDict[Sigmavec] = 0

    VisitedDict = Dict()
    counter = 0
    while counter < factorial(m - 1)
        counter += 1
        if mod(counter, 1000) == 0
            println("ronde ", counter)
        end
        #determine the current vertex
        currentvertex = nothing
        CurrentDistance = 1000000
        for (v, distance) in DistanceDict
            if distance < CurrentDistance && !haskey(VisitedDict, v)
                CurrentDistance = distance
                currentvertex = v
            end
        end

        VisitedDict[currentvertex] = true

        #iterate over neighbors of current vertex
        Nbrs = []
        #create neighbors of current vertex
        for i1 = 1:m
            i2 = i1 < m ? i1 + 1 : 1
            #make sure that i1=1 always
            if i2 == 1
                i2 = i1
                i1 = 1
            end
            NewNeighbor = deepcopy(currentvertex)
            NewNeighbor[[i1, i2]] = NewNeighbor[[i2, i1]]
            if i1 == 1
                NewNeighbor = vcat(NewNeighbor[i2:end], NewNeighbor[1:(i2-1)])
            end
            push!(Nbrs, NewNeighbor)
        end

        #Now iterate over the neighbors of current vertex
        for adjVertex in Nbrs
            if haskey(DistanceDict, adjVertex)
                if DistanceDict[adjVertex] > CurrentDistance + 1
                    DistanceDict[adjVertex] = CurrentDistance + 1
                end
            else
                DistanceDict[adjVertex] = CurrentDistance + 1
            end
        end
    end

    return DistanceDict
end




# attempt at faster implementation via the paper "Finding minimum-length generator sequences" by Mark R. Jerrum (Theoretical Computer Science 36 (1985))
# reference was suggested by Frank de Meijer 
function CsigmatauNew(m)

    @time Vertices = unique(labelCanonical(pushfirst!(vertex, 1)) for vertex in permutations(Array((2:m))))

    DistanceDict = Dict()  
    Xvec=zeros(Int,m)
    Ubound = m*m
    @time for vertex in Vertices
        minVal = Ubound 
        for shift in 0:(m-1) 
            vertexTemp = circshift(vertex,shift)
            for j=1:m 
                Xvec[j]=vertexTemp[j]-j
            end 
            while maximum(Xvec)-minimum(Xvec)>m 
                i=argmax(Xvec); j=argmin(Xvec); 
                Xvec[i] -= m;
                Xvec[j] += m; 
            end
            toCompare = IntersectionNumber(Xvec)
            toCompare < minVal ? minVal=toCompare : 0;
        end
        DistanceDict[vertex]=minVal 
    end
    return DistanceDict
end


function IntersectionNumber(xVector)
    m=size(xVector,1);
    Ival=0;
    for i=1:m, j=i+1:m 
        r=i-j;
        s=(i+xVector[i])-(j+xVector[j])
        if r <=s ## number of integers in [r,s] divisible by m
            Ival+=floor(Int,s/m)-ceil(Int,r/m) +1
        else #number of integers in [s,r] divisible by m
            Ival+= floor(Int,r/m)-ceil(Int,s/m)+1
        end
    end
    return Ival; 
end

#@time  d=Csigmatau(6)

#save("Crossing6.jld2", "CR6", d)
# @time Crossing9=load("Crossing9.jld2")["CR9"]

#println("full distances are generated for n=6")

# Attempt at a faster implementation. Not correct yet.
function Csigmatau2(m)

    DistanceDict = Dict()

    function getDist(pIn)#, isShifted = false)
        # first run: rotate cyclicly

        # Only gets called once, and sets the dictionary entry
        return get!(DistanceDict, pIn) do
            # @show pIn

            # ps = isShifted ? [pIn] : [circshift(pIn, i) for i = 1:m]
            ps = [pIn]
            #println(pIn)

            tmpRes = []

            for p in ps

                i = 1
                while findfirst(x -> x == i, p) == i && i < m
                    i += 1
                end

                if i == m
                    return 0
                end

                pos = findfirst(x -> x == i, p)

                forwardDist = pos - i
                backwardDist = m - pos + i - 1

                minDist = min(forwardDist, backwardDist)

                ind = [j >= i && j != pos for j = 1:m]

                newP = vcat(1:i, p[ind])
                # @show pos
                # @show p
                # @show i
                # @show forwardDist
                # @show backwardDist
                # @show newP

                push!(tmpRes, minDist + getDist(newP))#, !isShifted))
            end

            # @show tmpRes

            return minimum(tmpRes)
        end
    end

    # @show getDist([1,3,2,4])

    # first pass
    for s in permutations(1:m)
        # p = vcat([1],s)
        p = s

        if !haskey(DistanceDict, p)
            getDist(p)
        end
    end

    # second pass: circ shifts
    # not enough: shifting 1 later on may be better: 126345
    tmpDict = Dict()
    for s in permutations(2:m)
        p = vcat([1], s)

        tmpDict[p] = minimum([DistanceDict[circshift(p, i)] for i = 1:m])
    end

    # third pass: inversion
    resDict = Dict()
    for s in permutations(2:m)
        p = vcat([1], s)
        pInv = vcat([1], reverse(s))

        if tmpDict[p] != tmpDict[pInv]
            #display("AA")
        end

        resDict[p] = min(tmpDict[p], tmpDict[pInv])
    end

    return (tmpDict, resDict)
    return tmpDict
    return resDict
    # return DistanceDict
end


##

# k = 8
# @time A = Csigmatau(k)
# @time B = CsigmatauOld(k)

# println("B.count ",length(B))
# @show A== B

#  for k in keys(A)
#      if A[k] != B[k]
#          @show k
#          @show A[k]
#          @show B[k]
#      end
#  end