
using LinearAlgebra
include("CRRepresentativeSet.jl")

#Daniel's function to find a spanning set of columns fast
function findColSpanningSetV3(A)
    F = qr(A, Val(true))
    r = sum(abs.(diag(F.R)) .> 0.00000001)
    if r == 0
        return []
    end
    return sort(F.p[1:r])
end


function decomposeModuleBasicLinearAlgebra(m)
    #first we are going to consider the S_m-action on \C^{Z_m}
    blocksObject = generateBlocksObject(m)

    Sigmavec = Array((1:m))
    Vertices = collect(permutations(Array((2:m))))
    IndexDict = Dict()
    number = 1 
    for vertex in Vertices
        pushfirst!(vertex, 1)
        IndexDict[vertex]=number 
        number+=1
    end


    Blocks=[]
    test=0 
    @time for (key,val) in blocksObject
        A=zeros(factorial(m-1),size(val,1))
        Set= generateRepresentativeSetLambda10(val)
        index=1;
        for vector in Set
            v = zeros(factorial(m-1))
            for (cycle,value) in vector 
                v[IndexDict[cycle]]=value 
            end
            A[:,index]=v
            index+=1
        end
        SpanSet = findColSpanningSetV3(A)
        #println(SpanSet)
        mi= size(SpanSet,1)
        test+=mi^2
        if !isempty(SpanSet)
            push!(Blocks,(key,val[SpanSet]))
        end
    end

    #now we consider the additional S_2 action 
    println("------------------ NOW S2-part-------")
    BlocksAfterS2=[]
    test=0
    println("begin S2-part") 
    @time for (key,val) in Blocks
        A=zeros(factorial(m-1),size(val,1))
        Set= generateRepresentativeSetLambda10(val)
        ##first do +1 part	
        index=1;
        for vector in Set
            v = zeros(factorial(m-1))
            for (cycle,value) in vector
                v[IndexDict[cycle]]+=value
            cyclenew = deepcopy(cycle)
                cyclenew[2:end]=cyclenew[end:-1:2];
                v[IndexDict[cyclenew]]+=value
            end
            A[:,index]=v
            index+=1
        end
        SpanSet = findColSpanningSetV3(A)

        mi= size(SpanSet,1)
        test+=mi^2
        if !isempty(SpanSet)
            push!(BlocksAfterS2,((key,1),val[SpanSet]))
        end

        A=zeros(factorial(m-1),size(val,1))
        ##now do -1 part	
        index=1;
        for vector in Set
            v = zeros(factorial(m-1))
            for (cycle,value) in vector
                v[IndexDict[cycle]]+=value
                cyclenew = deepcopy(cycle)
                cyclenew[2:end]=cyclenew[end:-1:2];
                v[IndexDict[cyclenew]]-=value
            end
            A[:,index]=v
            index+=1
        end
        SpanSet = findColSpanningSetV3(A)
        mi= size(SpanSet,1)
        test+=mi^2
        if !isempty(SpanSet)
            push!(BlocksAfterS2,((key,-1),val[SpanSet]))
        end
    end


    mapLambdaToBlocksElement = Dict()
    for block in BlocksAfterS2
        mapLambdaToBlocksElement[block[1]] = block[2]
    end

    return mapLambdaToBlocksElement
end