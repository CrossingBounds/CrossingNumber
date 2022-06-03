using LinearAlgebra, AbstractAlgebra
using Combinatorics
using Printf
using JLD2, FileIO

setprecision(128)

#generate the young shapes of a tableau of k boxes and height at most r
function ShapeAtMost(k,r)
    startset = collect(partitions(k,1))

    for j=2:r
        startset = append!(startset,collect(partitions(k,j)))
    end

    return startset
end

function IsSemiStandard(Y)
    issemi=true;

    rowpartition = Y.part;
    colpartition = conj(Y).part;

    #check rows
    for i=1:size(rowpartition,1)
        for j=1:(rowpartition[i]-1);
            if (Y[i,j]> Y[i,j+1])
                return false
            end
        end
    end

    #check cols
    for i=1:size(colpartition,1)
        for j=1:(colpartition[i]-1);
            if (Y[j,i]>= Y[j+1,i])
                return false
            end
        end
    end

    return issemi;
end


function AllCandidateVectors(n)
    Candidates = permutations(collect(1:n));
    return Candidates
end


function ArrayCombinations(A,B)
    ret=[]
    for i=1:size(A,1)
        for j=1:size(B,1)
            ret=push!(ret,cat(A[i],B[j];dims=1))
        end
    end

    return ret;
end


#A is an array with signs, B is just an array of permutations
function ArrayCombinationsWithSigns(A,B)
    ret=[]
    for i=1:size(A,1)
        for j=1:size(B,1)
            ret=push!(ret,(cat(A[i][1],B[j];dims=1),A[i][2]*levicivita(B[j])))
        end
    end

    return ret;
end


#Returns an array with all young tableaux row equivalent to the given young tableau
function AllRowEquivalentTableaux(Y)
    RowEquivalentTableaux=[];
    rowpartition = Y.part;

    #begin with first row
    T=collect(multiset_permutations(Y[1,1:rowpartition[1]],rowpartition[1]));
    #add rows
    for i=2:size(rowpartition,1)
        T=ArrayCombinations(T, collect(multiset_permutations(Y[i,1:rowpartition[i]],rowpartition[i])))
    end

    for i=1:size(T,1)
        fill!(Y,T[i])
        RowEquivalentTableaux = push!(RowEquivalentTableaux,deepcopy(Y))
    end

    return RowEquivalentTableaux;
end


#Returns an array consisting of tuples (P*Y, sign P), where P runs over the permutations in the column stabilizer of Y
function AllColumnSignTableaux(Y)
    ColumnSignTableaux=[];
    rowpartition = Y.part;
    Yprime=conj(Y)
    colpartition = Yprime.part;

    #begin with first column, save permutation + sign
    Tstart=collect(permutations(collect(1:colpartition[1])));
    T=[];
    for i=1:size(Tstart,1)
        T=push!(T,(Tstart[i], levicivita(Tstart[i])))
    end

    #add columns, save permutation + sign
    for i=2:size(colpartition,1)
        T=ArrayCombinationsWithSigns(T, collect(permutations(collect(1:colpartition[i]))))
    end

    #make the new tableaux
    for i=1:size(T,1)
        start=1;
        Yprimecopy = deepcopy(Yprime)
        newfilling=[];
        for j=1:rowpartition[1]
            newfilling=append!(newfilling,Yprime[j,T[i][1][start:start+colpartition[j]-1]])
            start+=colpartition[j]
        end
        fill!(Yprimecopy,convert(Array{Int},newfilling))
        Ydash = conj(Yprimecopy)

        ColumnSignTableaux = push!(ColumnSignTableaux,(Ydash,T[i][2]))
    end

    return ColumnSignTableaux;
end



function generateBlocksObject(n)
    MapLambdaToBlocksElement = Dict();
    Lambdas = ShapeAtMost(n,n);
    blockSizes=[];
    for lambda in Lambdas
        #maxreprelementsize=0;
        Ystart = YoungTableau(lambda)
        GoodTableaux = []
        blockSize = 0;

        r=n;
        candidates = AllCandidateVectors(n);
        for candidate in candidates
            fill!(Ystart,candidate)
            if (IsSemiStandard(Ystart))
                #push a deepcopy so that we push the correct filling and do not change it afterwards.
                GoodTableaux = push!(GoodTableaux,deepcopy(Ystart))
            end
        end
        blockSize = size(GoodTableaux,1)

        if blockSize >0
            blockSizes=push!(blockSizes,blockSize);
            MapLambdaToBlocksElement[lambda] = GoodTableaux; 
        end

    end

    return MapLambdaToBlocksElement
end


function AllCandidateVectorsS2(n)
    test = permutations(collect(2:n-1));
    Candidates=[];
    for i in test
        Candidates =push!(Candidates,[[ones(Int,2,1); i]...]);
    end
    return Candidates
end

function generateBlocksObjectS2(n)
    MapLambdaToBlocksElement = Dict();
    Lambdas = ShapeAtMost(n,n);
    blockSizes=[];
    for lambda in Lambdas
        #maxreprelementsize=0;
        Ystart = YoungTableau(lambda)
        GoodTableaux = []
        blockSize = 0;

        r=n;
        candidates = AllCandidateVectorsS2(n);
        for candidate in candidates
            fill!(Ystart,candidate)
            if (IsSemiStandard(Ystart))
                #push a deepcopy so that we push the correct filling and do not change it afterwards.
                GoodTableaux = push!(GoodTableaux,deepcopy(Ystart))
            end
        end
        blockSize = size(GoodTableaux,1)

        if blockSize >0
            blockSizes=push!(blockSizes,blockSize);
            println("newblock")
            println(lambda)
            println(blockSize, " with square: ", blockSize*blockSize)
            println()
            MapLambdaToBlocksElement[lambda] = GoodTableaux; 
        end

    end

    return MapLambdaToBlocksElement
end



#generate representative set for the block indiced by lambda
function generateRepresentativeSetLambda(n,LambdaSSYTableaux)
    GoodTableaux = LambdaSSYTableaux
    blockSize=size(GoodTableaux,1)
    ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

    for rowindex = 1:blockSize
        sigma = GoodTableaux[rowindex];
        WordsWithSigns=[]
        RowTableaux = AllRowEquivalentTableaux(sigma)
        for rowtab in RowTableaux
            ColTableaux = AllColumnSignTableaux(rowtab)
            for coltab in ColTableaux
                FillVector = coltab[1].fill;
                #we combine the fillvector and the partition into an (n+1)-cycle: (1 + the permutation);
                Word = zeros(Int,n+1);
                Word[1]=1; 
                for symbolminone = 2:n+1
                    position = findall(x -> x.==symbolminone-1, FillVector)[1]
                    Word[symbolminone]=position+1;
                end
                Sign = coltab[2];
                WordsWithSigns=push!(WordsWithSigns, (Word, Sign))
            end
        end

        ReprArrayLambda = push!(ReprArrayLambda,WordsWithSigns)
    end

    return  ReprArrayLambda 
end

#generate column-representative set (without using the column-stabilizer) for the block indiced by lambda with corresponding array of SSYTs GoodTableaux
function generateColumnRepresentativeSetLambda(n,LambdaSSYTableaux)
    GoodTableaux = LambdaSSYTableaux
    blockSize=size(GoodTableaux,1)
    ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

    for rowindex = 1:blockSize
        sigma = GoodTableaux[rowindex];
        WordsWithSigns=[]
        RowTableaux = AllRowEquivalentTableaux(sigma)
        for rowtab in RowTableaux
            FillVector = rowtab.fill
            #we combine the fillvector and the partition into an (n+1)-cycle: (1 + the permutation);
            Word = zeros(Int,n+1);
            Word[1]=1; 
            for symbolminone = 2:n+1
                position = findall(x -> x.==symbolminone-1, FillVector)[1]
                Word[symbolminone]=position+1;
            end
            Sign = 1;
            WordsWithSigns=push!(WordsWithSigns, (Word, Sign))
        end

        ReprArrayLambda = push!(ReprArrayLambda,WordsWithSigns)
    end

    return  ReprArrayLambda 
end

# Word = [3,2,4,1]
# NewWord = deepcopy(Word)
# for i=1:length(Word)
#     NewWord[i]= findfirst(x-> x==i,Word)
# end
# println(NewWord)



function generateRepresentativeSetLambda10(GoodTableaux) 

    #BlocksObject= load("CrossingM7SymBasis.jld2")["basis"]
    #Output = Dict()
   # for (partition, GoodTableaux) in BlocksObject
        blockSize=size(GoodTableaux,1)
        ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

        for rowindex = 1:blockSize
            sigma = GoodTableaux[rowindex];
            WordsWithSigns=Dict{Array{Int8,1},Int64}();
            WordsWithSignsArray=Tuple{Array{Int8,1},Int64}[];
            RowTableaux = AllRowEquivalentTableaux(sigma)
            for rowtab in RowTableaux
                ColTableaux = AllColumnSignTableaux(rowtab)
                for coltab in ColTableaux
                    FillVector = coltab[1].fill;
                    #we combine the fillvector and the partition into an (n+1)-cycle: (1 + the permutation);
                    Word = [findfirst(x->x==i, FillVector) for i = 1:length(FillVector)]
                    indexwith1 = findfirst(x-> x==1,Word)
                    
                    Word = vcat(Word[indexwith1:end],Word[1:(indexwith1-1)])
                    
                    Sign = coltab[2];
                    
                    if haskey(WordsWithSigns,Word)
                        WordsWithSigns[Word]+= Sign
                    else
                        WordsWithSigns[Word]= Sign
                    end
                end
            end
            for (word,sign) in WordsWithSigns
                if sign !=0
                    push!(WordsWithSignsArray,(word,sign))
                end
            end

            ReprArrayLambda = push!(ReprArrayLambda,WordsWithSignsArray)
        end
        #Output[partition] = ReprArrayLambda
    #end
    #display(Output)
    return ReprArrayLambda
end
#generateRepresentativeSetLambda10() 

function generateColumnRepresentativeSetLambda10(GoodTableaux)  

    # BlocksObject= load("CrossingM7SymBasis.jld2")["basis"]
    # Output = Dict()
    # for (partition, GoodTableaux) in BlocksObject
        blockSize=size(GoodTableaux,1)
        ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

        for rowindex = 1:blockSize
            sigma = GoodTableaux[rowindex];
            WordsWithSigns=Dict{Array{Int8,1},Int64}();
            WordsWithSignsArray=Tuple{Array{Int8,1},Int64}[];
            RowTableaux = AllRowEquivalentTableaux(sigma)
            for rowtab in RowTableaux
                #ColTableaux = AllColumnSignTableaux(rowtab)
                #for coltab in ColTableaux
                    FillVector = rowtab.fill;
                    Word = [findfirst(x->x==i, FillVector) for i = 1:length(FillVector)]
                    
                    indexwith1 = findfirst(x-> x==1,Word)
                    Word = vcat(Word[indexwith1:end],Word[1:(indexwith1-1)])
                    
                    Sign = 1;

                    


                    if haskey(WordsWithSigns,Word)
                        WordsWithSigns[Word]+= Sign
                    else
                        WordsWithSigns[Word]= Sign
                    end
                #end
            end
            for (word,sign) in WordsWithSigns
                if sign !=0
                    push!(WordsWithSignsArray,(word,sign))
                end
            end

            ReprArrayLambda = push!(ReprArrayLambda,WordsWithSignsArray)
        end
        #Output[partition] = ReprArrayLambda
    #end
    return ReprArrayLambda#Output
end




function generateRepresentativeSetLambdaFULL(GoodTableaux, blocksign) 

    #BlocksObject= load("CrossingM7SymBasis.jld2")["basis"]
    #Output = Dict()
   # for (partition, GoodTableaux) in BlocksObject
        blockSize=size(GoodTableaux,1)
        ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

        for rowindex = 1:blockSize
            sigma = GoodTableaux[rowindex];
            WordsWithSigns=Dict{Array{Int8,1},Int64}();
            WordsWithSignsArray=Tuple{Array{Int8,1},Int64}[];
            RowTableaux = AllRowEquivalentTableaux(sigma)
            for rowtab in RowTableaux
                ColTableaux = AllColumnSignTableaux(rowtab)
                for coltab in ColTableaux
                    FillVector = coltab[1].fill;
                    #we combine the fillvector and the partition into an (n+1)-cycle: (1 + the permutation);
                    Word = [findfirst(x->x==i, FillVector) for i = 1:length(FillVector)]
                    indexwith1 = findfirst(x-> x==1,Word)
                    
                    Word = vcat(Word[indexwith1:end],Word[1:(indexwith1-1)])
                    
                    Word2 = deepcopy(Word)
                    Word2[2:end] = Word2[end:-1:2]

                    Sign = coltab[2];
                    if haskey(WordsWithSigns,Word)
                        WordsWithSigns[Word]+= Sign
                    else
                        WordsWithSigns[Word]= Sign
                    end

                    if blocksign == 1
                        if haskey(WordsWithSigns,Word2)
                            WordsWithSigns[Word2]+= Sign
                        else
                            WordsWithSigns[Word2]= Sign
                        end
                    end
                    if blocksign == -1
                        if haskey(WordsWithSigns,Word2)
                            WordsWithSigns[Word2]-= Sign
                        else
                            WordsWithSigns[Word2]= -Sign
                        end
                    end
                end
            end
            for (word,sign) in WordsWithSigns
                if sign !=0
                    push!(WordsWithSignsArray,(word,sign))
                end
            end

            ReprArrayLambda = push!(ReprArrayLambda,WordsWithSignsArray)
        end
        #Output[partition] = ReprArrayLambda
    #end
    #display(Output)
    return ReprArrayLambda
end
#generateRepresentativeSetLambda10() 

function generateColumnRepresentativeSetLambdaFULL(GoodTableaux, blocksign)  

    # BlocksObject= load("CrossingM7SymBasis.jld2")["basis"]
    # Output = Dict()
    # for (partition, GoodTableaux) in BlocksObject
        blockSize=size(GoodTableaux,1)
        ReprArrayLambda=Array{Tuple{Array{Int8,1},Int64}}[];

        for rowindex = 1:blockSize
            sigma = GoodTableaux[rowindex];
            WordsWithSigns=Dict{Array{Int8,1},Int64}();
            WordsWithSignsArray=Tuple{Array{Int8,1},Int64}[];
            RowTableaux = AllRowEquivalentTableaux(sigma)
            for rowtab in RowTableaux
                #ColTableaux = AllColumnSignTableaux(rowtab)
                #for coltab in ColTableaux
                    FillVector = rowtab.fill;
                    Word = [findfirst(x->x==i, FillVector) for i = 1:length(FillVector)]
                    
                    indexwith1 = findfirst(x-> x==1,Word)
                    Word = vcat(Word[indexwith1:end],Word[1:(indexwith1-1)])
                    
                    Sign = 1;

                    
                    Word2 = deepcopy(Word)
                    Word2[2:end] = Word2[end:-1:2]

                    if haskey(WordsWithSigns,Word)
                        WordsWithSigns[Word]+= Sign
                    else
                        WordsWithSigns[Word]= Sign
                    end

                    if blocksign == 1
                        if haskey(WordsWithSigns,Word2)
                            WordsWithSigns[Word2]+= Sign
                        else
                            WordsWithSigns[Word2]= Sign
                        end
                    end
                    if blocksign == -1
                        if haskey(WordsWithSigns,Word2)
                            WordsWithSigns[Word2]-= Sign
                        else
                            WordsWithSigns[Word2]= -Sign
                        end
                    end
                #end
            end
            for (word,sign) in WordsWithSigns
                if sign !=0
                    push!(WordsWithSignsArray,(word,sign))
                end
            end

            ReprArrayLambda = push!(ReprArrayLambda,WordsWithSignsArray)
        end
        #Output[partition] = ReprArrayLambda
    #end
    return ReprArrayLambda#Output
end





















#display(generateColumnRepresentativeSetLambda10())

function generateRepresentativeSetLambdaTestNew(n)#,LambdaSSYTableaux)
    BlocksObject= generateBlocksObject(n) #LambdaSSYTableaux
    Output = Dict()
    for (partition, GoodTableaux) in BlocksObject
        blockSize=size(GoodTableaux,1)
        ReprArrayLambda=[];

        for rowindex = 1:blockSize
            sigma = GoodTableaux[rowindex];
            WordsWithSigns=Dict();
            RowTableaux = AllRowEquivalentTableaux(sigma)
            for rowtab in RowTableaux
                ColTableaux = AllColumnSignTableaux(rowtab)
                for coltab in ColTableaux
                    FillVector = conj(coltab[1]).fill;
                    #we combine the fillvector and the partition into an (n+1)-cycle: (1 + the permutation);
                    Word = FillVector;

                    NewWord = deepcopy(Word)
                    for i=1:length(Word)
                        NewWord[i]= findfirst(x-> x==i,Word)
                    end
                    Word = NewWord
                    indexwith1 = findfirst(x-> x==1,Word)
                    Word = vcat(Word[indexwith1:end],Word[1:(indexwith1-1)])
                    
                    Sign = coltab[2];
                    
                    if haskey(WordsWithSigns,Word)
                        WordsWithSigns[Word]+= Sign
                    else
                        WordsWithSigns[Word]= Sign
                    end
                end
            end

            ReprArrayLambda = push!(ReprArrayLambda,WordsWithSigns)
        end
        Output[partition] = ReprArrayLambda
    end
    return Output
end
#display(generateRepresentativeSetLambdaTestNew(4))

#G = PermutationGroup(6)
#chi = character(Partition([2,2,2]))
#println(chi(perm"(1,2)(3,4)(5,6)"))