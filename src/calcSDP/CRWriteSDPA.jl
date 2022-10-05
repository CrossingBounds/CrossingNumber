using LinearAlgebra
using Combinatorics
using JLD2, FileIO
using Printf

setprecision(128)
include("CRRepresentativeSet.jl")
include("EquivalenceCR.jl")
include("Cmatrixsigmatau.jl")
include("BlockDiagonalize.jl")
#this is a test-comment to see if I can commit to GIT


#WRITE SEPARATE BLOCKS FOR ALPHA_M, CROSSING NUMBER PROBLEM. WE USE THE SYMMETRY GROUP S_{M-1} ACTING ON 
#Sigma_2,...,Sigma_{M} WHERE SIGMA IS AN M-CYCLE (seen as a row vector of length M)
#Where we assume that Sigma_1=1, i.e., {Sigma_2,...,Sigma_{M}} = {2,...,M}

function CrossingBlockDiagFULL(m)
    println("#### We consider alpha_m with m = ", m)
    println("### Using FULL reduction...")
    println("### Generating orbit numbers...")

    Numbers, maxnumber = EquivalenceClasses(m);

    ##### In order to speed-up lookups, we are going to store the numbers in an array
    ##### THIS COSTS A LOT OF MEMORY, SINCE (m-1)^(m-1) values are stored
    ##### BUT SIGNIFICANT SPEED INCREASE IN LOOKUPS
    # dimensionsNumbersArray= tuple([(m-1) for i=1:(m-1)]...)
    # #for m<=10, no numbers larger than 2^16 occur. So UInt16 is enough. For m=11, need UInt32. 
    # NumbersArray=zeros(UInt16, dimensionsNumbersArray);
    # for (cycle, number) in Numbers
    #     NumbersArray[[cycle[i]-1 for i =2:m]...]=number;
    # end
    NumbersDict = Dict(Numbers)

   # println(NumbersArray[3,2,1,4,5,6])
    #display(NumbersArray)
    ##### end of array creation



    println("largest orbit number (d_reduced) is equal to ", maxnumber)

    println("##### Generating Blocks...")
   # @time MapLambdaToBlocksElement =load(string("data/CrossingM",m,"FullSymBasis.jld2"))["basis"]
    
    @time MapLambdaToBlocksElement = decomposeModule(m)
    #display(MapLambdaToBlocksElement)
    #blockSizes = [factorial(m-1)]

    #println("block sizes: ",blockSizes)



    ##An auxiliary function which computes inner products. Need access to the dictionary Numbers, so we put it here.
    function ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
        EntryArray = zeros(Int,maxnumber)
        #println("inner product inner computation")
        @inbounds for cycle1WithCoef in ColReprSetEntry
            cycle1=cycle1WithCoef[1]

            #determine a dictionary corresponding with cycle1, (we are going to renumber cycle 1 to (1,...,m))
            #so that we can quickly renumber cycle2 to determine its orbit. 
            Cycle1Dict = deepcopy(cycle1)
            #Cycle1Dict[cycle1[1]] = 1 #just for reference, we don't use this value
            @inbounds for i=1:m 
                Cycle1Dict[cycle1[i]] = i
            end
            #Cycle1Dict=Cycle1Dict[2:m]
            #println(Cycle1Dict)
            #newc2 = deepcopy(cycle1); #can initialize it here, we will update it each loop. 
            @inbounds for cycle2WithCoef in RowReprSetEntry
                #renumber cycle2 to determine its orbit. 
                #coefficient of cycle 1 is always 1, so don't take it into account.
                #println(cycle2WithCoef[1])
                @inbounds Cycle2renumberedvector = Cycle1Dict[cycle2WithCoef[1]]

              #   println(Cycle2renumberedvector)
               # println(NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1])
              #  println(Cycle2renumberedvector)
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2:m] .- 1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                @inbounds EntryArray[NumbersDict[Cycle2renumberedvector]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[i]-1 for i = 2:m]] += cycle2WithCoef[2]*cycle1WithCoef[2];
            end
        end
        EntryDict=Dict()
        @inbounds for i=1:maxnumber
            if EntryArray[i]!=0
                EntryDict[i]=EntryArray[i]
            end
        end
        #println("inner product inner computation done")
        #println(EntryDict)
        return EntryDict;

    end

    for (partition, CorrespondingTableaux) in MapLambdaToBlocksElement
        println(size(CorrespondingTableaux,1))
    end
    println("##### Computing Block Diagonalization...")
    ReducedBlocks=[];
    blockSizes=[]
    @time for (partitionwithsign, CorrespondingTableaux) in MapLambdaToBlocksElement
        if partitionwithsign[1][1] != m  #we may disregard the partition [m], since it corresponds with the constraint that the sum of all entries in the matrix is nonnegative. We already have a linear equality <X,J>=1, so this constraint is superfluous. 
            blockSize =size(CorrespondingTableaux,1)
            push!(blockSizes,blockSize)
            println("##### Computing representative set for block indiced by ",partitionwithsign, " of size ",blockSize )
            RowReprSet = generateRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])
            ColReprSet = generateColumnRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])

            #println(RowReprSet)
            # println(ColReprSet)
            # println()
            blockSize = size(CorrespondingTableaux,1)

            NewBlock = Dict(); #we make a dictionary with NewBlock[[i,j]]= reduced entry. 
            for rowindex =1:blockSize
                println("rij ",rowindex)
                RowReprSetEntry = RowReprSet[rowindex]
                #println(RowReprSetEntry)
                for colindex=rowindex:blockSize
                    ColReprSetEntry = ColReprSet[colindex]
                    EntryDict = ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
                    NewBlock[[rowindex,colindex]] = EntryDict
                    # println("Entrydict")
                    # println(EntryDict)
                    # println()
                end
            end
            push!(ReducedBlocks,NewBlock)
       end
    end

    # return ReducedBlocks

    println("##### Writing the SDP...")
    ###start writing SDP file in SDPA-format (.dat-s)
    io = open(string("data/sdp_crossing_m$m","_S$m","xS2",".dat-s"), "w")

    # print number of variables and the block sizes
    # nVars = maxnumber 
    nVars = maxnumber; 
    println(io, nVars)

    #nBlocks = nVars (positivity constraints) + 2 (<X,J>=1 gives two inequality constraints) + number of blocks in the block diagonalization 
    nBlocks = nVars +2 + length(blockSizes);
    println(io, nBlocks)

    #print a line with the block sizes. First nVars (positivity constraints) + 2 (constraint <X,J>=1) times size 1
    for i=1:nVars+2
        print(io, 1," ")
    end
    #now the sizes of the blocks in the block diagonalization 
    for i=1:length(blockSizes)
        print(io, blockSizes[i]," ")
    end
    print(io,"\n")
    
    ### now we are going to construct the objective vector
    ### WE ARE GOING TO COMPUTE THE NECESSARY ENTRIES IN THE C-MATRIX,
    ### IN ORDER TO WRITE THE OBJECTIVE FUNCTION.
    cycles = collect(permutations(Array((2:m))))
    for cycle in cycles
        pushfirst!(cycle, 1)
    end

    distanceDict =Csigmatau(m)#load(string("Crossing",m,".jld2"))["CR$m"]##
    CmatrixDict = Dict() 
    #make distances from vertex to inverse of destination
    for (vertex, distance) in distanceDict 
        vertexnew = deepcopy(vertex)
        vertexnew[2:end]=vertexnew[end:-1:2];
        CmatrixDict[vertexnew] = distance
    end
    
    NumbersToCycleCoefficients = Dict()
    for i=1:maxnumber 
        NumbersToCycleCoefficients[i]=0;
    end
    for cycle in cycles 
        #we compute C((1,...,m), cycle). 
        currentOrbitNumber = Numbers[cycle]; 
        NumbersToCycleCoefficients[currentOrbitNumber] += CmatrixDict[cycle]
    end
    ## write down objective function 
    for i=1:maxnumber 
        print(io, NumbersToCycleCoefficients[i]," ")
    end
    print(io,"\n");

    #now print the nonnegativity constraints 
    for i=1:maxnumber 
        println(io,i," ",i," ",1, " ",1 , " ", 1)
    end
    currentblocknumber = maxnumber+1; 

    ### THE CONSTRAINT <X,J> = 1 is equivalent to (n-1)! * sum_{tau} X_{(1,..,n),tau} =1
    ### We will impose this constraint by first rewriting sum_{tau} X_{(1,..,n),tau} as a sum of y_i (orbit variables)
    NumbersToOrbitSizes = Dict()
    for i=1:maxnumber 
        NumbersToOrbitSizes[i]=0;
    end
    for cycle in cycles 
        currentOrbitNumber = Numbers[cycle]; 
        NumbersToOrbitSizes[currentOrbitNumber] += 1
    end
    # for i=1:maxnumber 
    #     print(NumbersToOrbitSizes[i]," ")
    # end

    # print the constraint sum alpha i y_i >=1, i.e. -1 + sum alpha_i y_i >=0 
    println(io,0," ",currentblocknumber," ",1, " ",1 , " ", 1)
    for i=1:maxnumber 
        println(io,i," ",currentblocknumber," ",1, " ",1 , " ", NumbersToOrbitSizes[i])  
    end
    currentblocknumber +=1
    # print the constraint sum alpha i y_i <=1, i.e. 1 - sum alpha_i y_i >=0 
    println(io,0," ",currentblocknumber," ",1, " ",1 , " ", -1)
    for i=1:maxnumber 
        println(io,i," ",currentblocknumber," ",1, " ",1 , " ", -NumbersToOrbitSizes[i])  
    end
    currentblocknumber +=1; 


    #Check block sizes 
    testsize = 0;
    for blockindex = 1:size(ReducedBlocks,1)
        blockSize = blockSizes[blockindex]
        testsize += blockSize*blockSize;
    end
    println("#### Sum of squares of block sizes: ", testsize)

    ### NOW WRITE BLOCK DIAGONALIZATION. 
    for blockindex = 1:size(ReducedBlocks,1)
        blockSize = blockSizes[blockindex]
        Block = ReducedBlocks[blockindex]

        for i=1:blockSize 
            for j=i:blockSize 
                EntryDict = Block[[i,j]]
                for varNumber = 1:nVars 
                    if haskey(EntryDict,varNumber) #&& EntryDict[varNumber] != 0
                        println(io,varNumber," ",currentblocknumber," ",i, " ",j , " ", EntryDict[varNumber])  
                    end
                end
            end
        end

        currentblocknumber+=1; 
    end

    ### HERE IS CODE WHICH WAS USED TO WRITE THE FULL MOMENT MATRIX WITHOUT BLOCK DIAGONALIZATION TO EXPERIMENT. 

    # for rowidx = 1:factorial(m-1)
    #     sigma = cycles[rowidx]; 
    # for colidx = rowidx:factorial(m-1)
    #     tau = cycles[colidx]  
    #     (newsigma, newtau) = RenumberPermutations(sigma,tau)
    #     println(io,Numbers[newtau]," ",currentblocknumber," ",rowidx, " ",colidx , " ", 1)  
    # end
    # end
    close(io)
    println("##### Writing finished.")
end





function CrossingBlockDiagOnly11block(m)
    println("#### We consider alpha_m with m = ", m)
    println("### Using FULL reduction...")
    println("### Generating orbit numbers...")

    Numbers, maxnumber = EquivalenceClasses(m);

    ##### In order to speed-up lookups, we are going to store the numbers in an array
    ##### THIS COSTS A LOT OF MEMORY, SINCE (m-1)^(m-1) values are stored
    ##### BUT SIGNIFICANT SPEED INCREASE IN LOOKUPS
    # dimensionsNumbersArray= tuple([(m-1) for i=1:(m-1)]...)
    # #for m<=10, no numbers larger than 2^16 occur. So UInt16 is enough. For m=11, need UInt32. 
    # NumbersArray=zeros(UInt16, dimensionsNumbersArray);
    # for (cycle, number) in Numbers
    #     NumbersArray[[cycle[i]-1 for i =2:m]...]=number;
    # end
    NumbersDict = Dict(Numbers)

   # println(NumbersArray[3,2,1,4,5,6])
    #display(NumbersArray)
    ##### end of array creation



    println("largest orbit number (d_reduced) is equal to ", maxnumber)

    println("##### Generating Blocks...")
   # @time MapLambdaToBlocksElement =load(string("data/CrossingM",m,"FullSymBasis.jld2"))["basis"]
    
    @time MapLambdaToBlocksElement = decomposeModule(m)

    for (key,value) in MapLambdaToBlocksElement
        if !(key[1][1]==m-2 && key[1][2]==1)
            delete!(MapLambdaToBlocksElement,key)
        end
    end
    #display(MapLambdaToBlocksElement)
    #blockSizes = [factorial(m-1)]

    #println("block sizes: ",blockSizes)



    ##An auxiliary function which computes inner products. Need access to the dictionary Numbers, so we put it here.
    function ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
        EntryArray = zeros(Int,maxnumber)
        #println("inner product inner computation")
        @inbounds for cycle1WithCoef in ColReprSetEntry
            cycle1=cycle1WithCoef[1]

            #determine a dictionary corresponding with cycle1, (we are going to renumber cycle 1 to (1,...,m))
            #so that we can quickly renumber cycle2 to determine its orbit. 
            Cycle1Dict = deepcopy(cycle1)
            #Cycle1Dict[cycle1[1]] = 1 #just for reference, we don't use this value
            @inbounds for i=1:m 
                Cycle1Dict[cycle1[i]] = i
            end
            #Cycle1Dict=Cycle1Dict[2:m]
            #println(Cycle1Dict)
            #newc2 = deepcopy(cycle1); #can initialize it here, we will update it each loop. 
            @inbounds for cycle2WithCoef in RowReprSetEntry
                #renumber cycle2 to determine its orbit. 
                #coefficient of cycle 1 is always 1, so don't take it into account.
                #println(cycle2WithCoef[1])
                @inbounds Cycle2renumberedvector = Cycle1Dict[cycle2WithCoef[1]]

              #   println(Cycle2renumberedvector)
               # println(NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1])
              #  println(Cycle2renumberedvector)
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2:m] .- 1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                @inbounds EntryArray[NumbersDict[Cycle2renumberedvector]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[i]-1 for i = 2:m]] += cycle2WithCoef[2]*cycle1WithCoef[2];
            end
        end
        EntryDict=Dict()
        @inbounds for i=1:maxnumber
            if EntryArray[i]!=0
                EntryDict[i]=EntryArray[i]
            end
        end
        #println("inner product inner computation done")
        #println(EntryDict)
        return EntryDict;

    end

    for (partition, CorrespondingTableaux) in MapLambdaToBlocksElement
        println(size(CorrespondingTableaux,1))
    end
    println("##### Computing Block Diagonalization...")
    ReducedBlocks=[];
    blockSizes=[]
    @time for (partitionwithsign, CorrespondingTableaux) in MapLambdaToBlocksElement
        if partitionwithsign[1][1] != m  #we may disregard the partition [m], since it corresponds with the constraint that the sum of all entries in the matrix is nonnegative. We already have a linear equality <X,J>=1, so this constraint is superfluous. 
            blockSize =size(CorrespondingTableaux,1)
            push!(blockSizes,blockSize)
            println("##### Computing representative set for block indiced by ",partitionwithsign, " of size ",blockSize )
            RowReprSet = generateRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])
            ColReprSet = generateColumnRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])

            #println(RowReprSet)
            # println(ColReprSet)
            # println()
            blockSize = size(CorrespondingTableaux,1)

            NewBlock = Dict(); #we make a dictionary with NewBlock[[i,j]]= reduced entry. 
            for rowindex =1:blockSize
                println("rij ",rowindex)
                RowReprSetEntry = RowReprSet[rowindex]
                #println(RowReprSetEntry)
                for colindex=rowindex:blockSize
                    ColReprSetEntry = ColReprSet[colindex]
                    EntryDict = ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
                    NewBlock[[rowindex,colindex]] = EntryDict
                    # println("Entrydict")
                    # println(EntryDict)
                    # println()
                end
            end
            push!(ReducedBlocks,NewBlock)
       end
    end

    # return ReducedBlocks

    println("##### Writing the SDP...")
    ###start writing SDP file in SDPA-format (.dat-s)
    io = open(string("data/sdp_crossingConjectureBlock_m$m","_S$m","xS2",".dat-s"), "w")
    # print number of variables and the block sizes
    # nVars = maxnumber 
    nVars = maxnumber; 
    println(io, nVars)

    #nBlocks = nVars (positivity constraints) + 2 (<X,J>=1 gives two inequality constraints) + number of blocks in the block diagonalization 
    nBlocks = nVars +2 + length(blockSizes);
    println(io, nBlocks)

    #print a line with the block sizes. First nVars (positivity constraints) + 2 (constraint <X,J>=1) times size 1
    for i=1:nVars+2
        print(io, 1," ")
    end
    #now the sizes of the blocks in the block diagonalization 
    for i=1:length(blockSizes)
        print(io, blockSizes[i]," ")
    end
    print(io,"\n")
    
    ### now we are going to construct the objective vector
    ### WE ARE GOING TO COMPUTE THE NECESSARY ENTRIES IN THE C-MATRIX,
    ### IN ORDER TO WRITE THE OBJECTIVE FUNCTION.
    cycles = collect(permutations(Array((2:m))))
    for cycle in cycles
        pushfirst!(cycle, 1)
    end

    distanceDict =Csigmatau(m)#load(string("Crossing",m,".jld2"))["CR$m"]##
    CmatrixDict = Dict() 
    #make distances from vertex to inverse of destination
    for (vertex, distance) in distanceDict 
        vertexnew = deepcopy(vertex)
        vertexnew[2:end]=vertexnew[end:-1:2];
        CmatrixDict[vertexnew] = distance
    end
    
    NumbersToCycleCoefficients = Dict()
    for i=1:maxnumber 
        NumbersToCycleCoefficients[i]=0;
    end
    for cycle in cycles 
        #we compute C((1,...,m), cycle). 
        currentOrbitNumber = Numbers[cycle]; 
        NumbersToCycleCoefficients[currentOrbitNumber] += CmatrixDict[cycle]
    end
    ## write down objective function 
    for i=1:maxnumber 
        print(io, NumbersToCycleCoefficients[i]," ")
    end
    print(io,"\n");

    #now print the nonnegativity constraints 
    for i=1:maxnumber 
        println(io,i," ",i," ",1, " ",1 , " ", 1)
    end
    currentblocknumber = maxnumber+1; 

    ### THE CONSTRAINT <X,J> = 1 is equivalent to (n-1)! * sum_{tau} X_{(1,..,n),tau} =1
    ### We will impose this constraint by first rewriting sum_{tau} X_{(1,..,n),tau} as a sum of y_i (orbit variables)
    NumbersToOrbitSizes = Dict()
    for i=1:maxnumber 
        NumbersToOrbitSizes[i]=0;
    end
    for cycle in cycles 
        currentOrbitNumber = Numbers[cycle]; 
        NumbersToOrbitSizes[currentOrbitNumber] += 1
    end
    # for i=1:maxnumber 
    #     print(NumbersToOrbitSizes[i]," ")
    # end

    # print the constraint sum alpha i y_i >=1, i.e. -1 + sum alpha_i y_i >=0 
    println(io,0," ",currentblocknumber," ",1, " ",1 , " ", 1)
    for i=1:maxnumber 
        println(io,i," ",currentblocknumber," ",1, " ",1 , " ", NumbersToOrbitSizes[i])  
    end
    currentblocknumber +=1
    # print the constraint sum alpha i y_i <=1, i.e. 1 - sum alpha_i y_i >=0 
    println(io,0," ",currentblocknumber," ",1, " ",1 , " ", -1)
    for i=1:maxnumber 
        println(io,i," ",currentblocknumber," ",1, " ",1 , " ", -NumbersToOrbitSizes[i])  
    end
    currentblocknumber +=1; 


    #Check block sizes 
    testsize = 0;
    for blockindex = 1:size(ReducedBlocks,1)
        blockSize = blockSizes[blockindex]
        testsize += blockSize*blockSize;
    end
    println("#### Sum of squares of block sizes: ", testsize)

    ### NOW WRITE BLOCK DIAGONALIZATION. 
    for blockindex = 1:size(ReducedBlocks,1)
        blockSize = blockSizes[blockindex]
        Block = ReducedBlocks[blockindex]

        for i=1:blockSize 
            for j=i:blockSize 
                EntryDict = Block[[i,j]]
                for varNumber = 1:nVars 
                    if haskey(EntryDict,varNumber) #&& EntryDict[varNumber] != 0
                        println(io,varNumber," ",currentblocknumber," ",i, " ",j , " ", EntryDict[varNumber])  
                    end
                end
            end
        end

        currentblocknumber+=1; 
    end
    close(io)
    println("##### Writing finished.")
end



# test = CrossingBlockDiagFULL(7)

#@time MapLambdaToBlocksElement =load(string("CrossingM",8,"SymBasis.jld2"))["basis"]

#println(MapLambdaToBlocksElement)

#@time CrossingBlockDiag(4)


# @time M=zeros( 9,9,9,9,9,9,9,9,9);
# println(M[1,1,1,1,1,1,1,1,1])
# println("M success")

# m=9
# test= tuple(vcat([1],[m for i=1:m])...)
# println(test)


## 

function CrossingBlockDiagOnly(m)
    println("#### We consider alpha_m with m = ", m)
    println("### Using FULL reduction...")
    println("### Generating orbit numbers...")

    Numbers, maxnumber = EquivalenceClasses(m);

    ##### In order to speed-up lookups, we are going to store the numbers in an array
    ##### THIS COSTS A LOT OF MEMORY, SINCE (m-1)^(m-1) values are stored
    ##### BUT SIGNIFICANT SPEED INCREASE IN LOOKUPS
    # dimensionsNumbersArray= tuple([(m-1) for i=1:(m-1)]...)
    # #for m<=10, no numbers larger than 2^16 occur. So UInt16 is enough. For m=11, need UInt32. 
    # NumbersArray=zeros(UInt16, dimensionsNumbersArray);
    # for (cycle, number) in Numbers
    #     NumbersArray[[cycle[i]-1 for i =2:m]...]=number;
    # end
    NumbersDict = Dict(Numbers)

   # println(NumbersArray[3,2,1,4,5,6])
    #display(NumbersArray)
    ##### end of array creation



    println("largest orbit number (d_reduced) is equal to ", maxnumber)

    println("##### Generating Blocks...")
    @time MapLambdaToBlocksElement = decomposeModule(9)# load("data/decompose9.jld2")["basis"]
    #display(MapLambdaToBlocksElement)
    #blockSizes = [factorial(m-1)]

    #println("block sizes: ",blockSizes)



    ##An auxiliary function which computes inner products. Need access to the dictionary Numbers, so we put it here.
    function ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
        EntryArray = zeros(Int,maxnumber)
        #println("inner product inner computation")
        @inbounds for cycle1WithCoef in ColReprSetEntry
            cycle1=cycle1WithCoef[1]

            #determine a dictionary corresponding with cycle1, (we are going to renumber cycle 1 to (1,...,m))
            #so that we can quickly renumber cycle2 to determine its orbit. 
            Cycle1Dict = deepcopy(cycle1)
            #Cycle1Dict[cycle1[1]] = 1 #just for reference, we don't use this value
            @inbounds for i=1:m 
                Cycle1Dict[cycle1[i]] = i
            end
            #Cycle1Dict=Cycle1Dict[2:m]
            #println(Cycle1Dict)
            #newc2 = deepcopy(cycle1); #can initialize it here, we will update it each loop. 
            @inbounds for cycle2WithCoef in RowReprSetEntry
                #renumber cycle2 to determine its orbit. 
                #coefficient of cycle 1 is always 1, so don't take it into account.
                #println(cycle2WithCoef[1])
                @inbounds Cycle2renumberedvector = Cycle1Dict[cycle2WithCoef[1]]

              #   println(Cycle2renumberedvector)
               # println(NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1])
              #  println(Cycle2renumberedvector)
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2]-1,Cycle2renumberedvector[3]-1,Cycle2renumberedvector[4]-1,Cycle2renumberedvector[5]-1,Cycle2renumberedvector[6]-1,Cycle2renumberedvector[7]-1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[2:m] .- 1]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                @inbounds EntryArray[NumbersDict[Cycle2renumberedvector]] += cycle2WithCoef[2]*cycle1WithCoef[2];
                
                # @inbounds EntryArray[NumbersArray[Cycle2renumberedvector[i]-1 for i = 2:m]] += cycle2WithCoef[2]*cycle1WithCoef[2];
            end
        end
        EntryDict=Dict()
        @inbounds for i=1:maxnumber
            if EntryArray[i]!=0
                EntryDict[i]=EntryArray[i]
            end
        end
        #println("inner product inner computation done")
        #println(EntryDict)
        return EntryDict;

    end

    for (partition, CorrespondingTableaux) in MapLambdaToBlocksElement
        println(size(CorrespondingTableaux,1))
    end
    println("##### Computing Block Diagonalization...")
    ReducedBlocks=[];
    blockSizes=[]
    @time for (partitionwithsign, CorrespondingTableaux) in MapLambdaToBlocksElement
        if true# partitionwithsign[1].part == [6,3] #partitionwithsign[1][1] != m  #we may disregard the partition [m], since it corresponds with the constraint that the sum of all entries in the matrix is nonnegative. We already have a linear equality <X,J>=1, so this constraint is superfluous. 
            blockSize =size(CorrespondingTableaux,1)
            push!(blockSizes,blockSize)
            println("##### Computing representative set for block indiced by ",partitionwithsign, " of size ",blockSize )
            RowReprSet = generateRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])
            ColReprSet = generateColumnRepresentativeSetLambdaFULL(CorrespondingTableaux, partitionwithsign[2])

            #println(RowReprSet)
            # println(ColReprSet)
            # println()
            blockSize = size(CorrespondingTableaux,1)

            NewBlock = Dict(); #we make a dictionary with NewBlock[[i,j]]= reduced entry. 
            for rowindex =1:blockSize
                println("rij ",rowindex)
                RowReprSetEntry = RowReprSet[rowindex]
                #println(RowReprSetEntry)
                for colindex=rowindex:blockSize
                    ColReprSetEntry = ColReprSet[colindex]
                    EntryDict = ComputeCrossingInnerProduct(RowReprSetEntry, ColReprSetEntry)
                    NewBlock[[rowindex,colindex]] = EntryDict
                    # println("Entrydict")
                    # println(EntryDict)
                    # println()
                end
            end
            push!(ReducedBlocks,NewBlock)
       end
    end

    return (NumbersDict, ReducedBlocks)
end
