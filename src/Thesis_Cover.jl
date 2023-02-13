include("calcSDP/BuildSDP.jl")


m=8

using LinearAlgebra, MosekTools
# using SDPAFamily, LinearAlgebra#,JLD2

mu = (Partition([m - 2, 1, 1]), -1)
dualBlocksRestricted = [mu]


# sdpData = load("data/sdpDataBeta$m.jld2")["sdpData"]

@time sdpData = calcSDP(m, dualBlocksRestricted);
# save("data/sdpDataBeta$m.jld2", "sdpData", sdpData)

@time (t, yD, cD) = buildDualSDP_GMP(sdpData, dualBlocksRestricted);

addedConstraints = Set()
addedConstraintsDict = Dict()
addedConstraintsInd = Set()
bigM = 1e+5
constraints = Any[t<=bigM]
constraints = Any[c[2] >= 0 for c in cD]


for m in keys(yD)
    push!(constraints, sum(yD[m]) >= -100)
end

it = 1
maxIt = 100

setprecision(256)

@time while true
    if it == maxIt
        @warn "Reached iteration limit, aborting..."
        break
    end
    @info("Iteration $(it)")
    it += 1

    # P = maximize(t, constraints...; numeric_type=BigFloat)
    # P = maximize(0, constraints...; numeric_type=BigFloat)
    # P = maximize(t, constraints...; numeric_type=Float64) 
    P = maximize(0, constraints...; numeric_type=Float64)

    # params = SDPAFamily.Params{:sdpa_gmp,BigFloat}(
    #     maxIteration=200,
    #     epsilonStar=1e-40,
    #     lambdaStar=1e5,
    #     omegaStar=2,
    #     lowerBound=-1e5,
    #     upperBound=1e5,
    #     betaStar=0.1,
    #     betaBar=0.2,
    #     gammaStar=0.9,
    #     epsilonDash=1e-40,
    #     precision=512,
    # )

    @info "starting solver"

    # @time solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=params),silent_solver=false)
    # @time solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=SDPAFamily.UNSTABLE_BUT_FAST),silent_solver=false)
    @time solve!(P, () -> Mosek.Optimizer())

    @show P.optval
    @show P.status

    global optVal = P.optval

    maxViolation = 1
    ind = 0
    lastY = Convex.evaluate(yD[mu])
    lastT = Convex.evaluate(t)

    
    @showprogress "checking constraints" for (i,v) in enumerate(sdpData.vars)
        if !(i in addedConstraintsInd)
            constr = 0
            if haskey(sdpData.coeff[mu],[i for i in v])
                constr = -dot(Symmetric(sdpData.coeff[mu][[i for i in v]]), lastY)
            end
            constr += - lastT*sdpData.orbitSizes[v] + sdpData.obj[v]

            if constr < maxViolation
                maxViolation = constr 
                ind = i
            end
        end
    end


    @show maxViolation
    if true#maxViolation >= 0.000001
        @info("SDP Solved!")
        break
    else
        @info("Adding constraint $ind")
        push!(addedConstraints, cD[ind])
        addedConstraintsDict[ind] = cD[ind][2] >= 0
        push!(constraints, cD[ind][2] >= 0)
        push!(addedConstraintsInd, ind)
    end
end

@show 8 * optVal / (m * (m - 1))

# opt = Symmetric(Convex.evaluate(yD[mu]))
cent = Symmetric(Convex.evaluate(yD[mu]))
##

d = load("preciseSolution.jld2")
opt = d["opt"]
cent = d["cent"]

##

using MosekTools, JuMP

m = Model(Mosek.Optimizer)
set_silent(m)

obj = @variable(m,obj)
t = @variable(m,t)

ycent = opt + 1.0*(cent-opt)

mS = size(opt,1)

dirY = opt - ycent
startY = @variable(m, startY[1:mS, 1:mS], Symmetric)

Y = startY + obj * dirY
@constraint(m, Y in PSDCone())

for (o, c) in sdpData.coeff[mu]
    @constraint(m, -dot(c, Y) - t*sdpData.orbitSizes[Tuple(o)] + sdpData.obj[Tuple(o)] >= 0)
end

@objective(m, Max, t)
# P = maximize(t, constr...; numeric_type=Float64)

# solve!(P, () -> Mosek.Optimizer(), silent_solver=true)



# @show [Convex.evaluate(c[2]) for c in cD]


function distanceTillFeasible(stY)

    # if minimum(eigen(Float64.(stY)).values) < 0
    #     return -100
    # end

    fix.(startY, Symmetric(Float64.(stY)), force=true)

    optimize!(m)

    # @show P.status
    
    # @show termination_status(m)

    if termination_status(m) in [JuMP.MathOptInterface.INFEASIBLE]#,JuMP.MathOptInterface.SLOW_PROGRESS]
        return -100
    end
 
    return value(t)
end
# function distanceTillFeasible(startT, startY, dirT, dirY)

#     if minimum(eigen(Float64.(startY)).values) < 0
#         return -100
#     end

#     obj = Variable()
#     t = Variable()

#     Y = startY + obj * dirY
#     constr = Any[Y in :SDP]

#     for (o, c) in sdpData.coeff[mu]
#         push!(constr, -dot(c, Y) - t*sdpData.orbitSizes[Tuple(o)] + sdpData.obj[Tuple(o)] >= 0)
#     end

#     P = maximize(t, constr...; numeric_type=Float64)

#     solve!(P, () -> Mosek.Optimizer(), silent_solver=true)

#     @show P.status

#     if P.status in [JuMP.MathOptInterface.INFEASIBLE]#,JuMP.MathOptInterface.SLOW_PROGRESS]
#         return -100
#     end

#     # @show [Convex.evaluate(c[2]) for c in cD]

#     return Convex.evaluate(t)
# end
# function distanceTillFeasible(startT, startY, dirT, dirY)

#     if minimum(eigen(Float64.(startY)).values) < 0
#         return -100
#     end

#     obj = Variable()

#     constr = Any[c[2] >= 0 for c in cD]
#     # push!(constr, t == startT + obj *dirT)
#     push!(constr, yD[mu] == startY + obj * dirY)
#     push!(constr, yD[mu] == yD[mu]')
    
#     # P = maximize(obj, constr...; numeric_type=BigFloat)
#     P = maximize(t, constr...; numeric_type=Float64)

#     params = SDPAFamily.Params{:sdpa_gmp,BigFloat}(
#         maxIteration=200,
#         epsilonStar=1e-40,
#         lambdaStar=1e5,
#         omegaStar=2,
#         lowerBound=-1e5,
#         upperBound=1e5,
#         betaStar=0.1,
#         betaBar=0.2,
#         gammaStar=0.9,
#         epsilonDash=1e-40,
#         precision=64,
#     )
#     # solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=params))
#     # solve!(P, () -> SDPAFamily.Optimizer(presolve=true, params=SDPAFamily.UNSTABLE_BUT_FAST))
#     solve!(P, () -> Mosek.Optimizer())

#     @show P.status

#     if P.status in [JuMP.MathOptInterface.INFEASIBLE]#,JuMP.MathOptInterface.SLOW_PROGRESS]
#         return -100
#     end

#     @show [Convex.evaluate(c[2]) for c in cD]

#     return Convex.evaluate(t)
# end

# tStart = 0
# tDir = optVal

# yStart = opt#[0 0; 0 0]
# yDir = [0 0; 0 0]#opt

# distanceTillFeasible(tStart, yStart, tDir, yDir)


##

function randPSD(n)
    res = zeros(n,n)
    for i = 1:n
        v = 2*rand(n) .- 1
        res += rand()*v*v'
    end
    return Symmetric(res)
end

mS = size(opt, 1)


counter = 1
while true
    # ycent = 0.1*randPSD(mS) + opt

    # lambda = 0.1
    # ycent = opt + 0.5*(cent-opt)
    eigen(Float64.(ycent)).values
    xDir = Symmetric(2 .*rand(mS, mS) .- 1)
    xDir -= (opt - ycent)*dot(xDir, (opt - ycent))/dot((opt - ycent), (opt - ycent))
    # xDir = 0.004*xDir ./ dot(xDir, xDir)

    # for m = 9
    # xDir = 0.00008*xDir ./ dot(xDir, xDir)
    # for m = 8
    xDir = 0.0003*xDir ./ dot(xDir, xDir)
    # for m = 7
    # xDir = 0.0008*xDir ./ dot(xDir, xDir)

    yDir = Symmetric(2 .*rand(mS, mS) .- 1)
    yDir -= (opt - ycent)*dot(yDir, (opt - ycent))/dot((opt - ycent), (opt - ycent))
    yDir -= xDir * dot(yDir, xDir) / dot(xDir, xDir)
    # yDir = 0.004*yDir ./ dot(yDir, yDir)

    # for m = 9
    # yDir = 0.00008*yDir ./ dot(yDir, yDir)
    # for m = 8
    yDir = 0.0003*yDir ./ dot(yDir, yDir)
    # for m = 7
    # yDir = 0.0008*yDir ./ dot(yDir, yDir)

    using JLD2
    # save("GoodDirsm9.jld2", "ycent", ycent, "xDir", xDir, "yDir", yDir)

    #
    res = 100

    # @time picture = [distanceTillFeasible(0, ycent + x*xDir + y*yDir, 1, opt - ycent) for x=-1:(1/res):1, y=-1:(1/res):1]
    @time picture = [distanceTillFeasible(ycent + x*xDir + y*yDir) for x=-1:(1/res):1, y=-1:(1/res):1]

    # distanceTillFeasible(0, Float64.(ycent), 1, Float64.(opt - ycent))
    # distanceTillFeasible(Float64.(ycent))

    # Old, boring variant
    # ycent = 0.1*randPSD(mS) + opt
    # ycent = opt + 0.5*(cent-opt)
    # xDir = Symmetric(2 .*rand(mS, mS) .- 1)
    # xDir -= ycent*dot(xDir, ycent)/dot(ycent, ycent)
    # xDir = 0.001*xDir ./ dot(xDir, xDir)

    # yDir = Symmetric(2 .*rand(mS, mS) .- 1)
    # yDir -= ycent*dot(yDir, ycent)/dot(ycent, ycent)
    # yDir = 0.001*yDir ./ dot(yDir, yDir)


    # res = 10

    # picture = [distanceTillFeasible(0, ycent + x*xDir + y*yDir, 1, zero(ycent)) for x=-1:(1/res):1, y=-1:(1/res):1]

    using Plots

    pic = deepcopy(picture)
    pic[pic .<= -100] .= minimum(pic[pic .> -100])

    heatmap(Float64.(pic))

    # distanceTillFeasible(0, ycent, 1, zero(ycent))

    #
    gr()
    p = plot(Float64.(pic), st=:surface)

    
    save("pic$(counter).png", p)
    
    surfaceDir = [[0.0,0.0,0.0] for i = 1:size(picture,1)-1,j=1:size(picture,1)-1]
    
    for (ix, x) in enumerate(-1:(1/res):1)
        for (iy, y) in enumerate(-1:(1/res):1)
            if ix < size(picture,1) && iy < size(picture,1)
                a = [x,y,picture[ix,iy]]
                b = [x+1/res,y,picture[ix+1,iy]]
                c = [x,y+1/res,picture[ix,iy+1]]
                surfaceDir[ix,iy] .= cross(b-a, c-a)
                surfaceDir[ix, iy] ./= sqrt(dot(surfaceDir[ix, iy],surfaceDir[ix, iy]))
                surfaceDir[ix, iy] .+= 1
                surfaceDir[ix, iy] ./= 2.1
            end
        end
    end
    
    pictureNoise = [perlin(s...) for s in surfaceDir]
    p = heatmap(pictureNoise)
    
    # save("picNoise$(counter).png", p)
    # save("pic$(counter).jld2", "ycent", ycent, "xDir", xDir, "yDir", yDir, "picture", picture, "pictureNoise", pictureNoise)
    
    counter += 1
end

##
# picNo = 39
# picNo = 42
picNo = 8
# picNo = 54
# picNo = 52
# picNo = 41
# picNo = 23
# picNo = 19



# for picNo in [39,42,8,19]

data = load("pic$(picNo).jld2")

xDir = data["xDir"] * 1.6
yDir = data["yDir"] * 1.6
ycent = data["yCenter"]

res = 1000

function compImage(yc, xD, yD, r)

    P = zeros(Float64, 2*res+1, 2*res+1)

    @showprogress for (i, x) in enumerate(-1:(1/r):1)
        for (j, y) in enumerate(-1:(1/r):1)
            P[i,j] = distanceTillFeasible(yc + x*xD + y*yD)
        end
    end

    # @time picture = [distanceTillFeasible(yc + x*xD + y*yD) for x=-1:(1/r):1, y=-1:(1/r):1]
    return P
end

@time picture = compImage(ycent, xDir, yDir, res)

save("backData$(picNo).jld2", "yCenter", ycent, "xDir", xDir, "yDir", yDir, "picture", picture)

# @time picture = [distanceTillFeasible(ycent + x*xDir + y*yDir) for x=-1:(1/res):1, y=-1:(1/res):1]

using Plots

pic = deepcopy(picture)
pic[pic .<= -100] .= minimum(pic[pic .> -100])

# heatmap(Float64.(pic))

# distanceTillFeasible(0, ycent, 1, zero(ycent))

#
gr()
# p = plot(Float64.(pic), st=:surface)


# save("pic$(counter).png", p)

surfaceDir = [[0.0,0.0,0.0] for i = 1:size(picture,1)-1,j=1:size(picture,1)-1]

for (ix, x) in enumerate(-1:(1/res):1)
    for (iy, y) in enumerate(-1:(1/res):1)
        if ix < size(picture,1) && iy < size(picture,1)
            a = [x,y,picture[ix,iy]]
            b = [x+1/res,y,picture[ix+1,iy]]
            c = [x,y+1/res,picture[ix,iy+1]]
            surfaceDir[ix,iy] .= cross(b-a, c-a)
            surfaceDir[ix, iy] ./= sqrt(dot(surfaceDir[ix, iy],surfaceDir[ix, iy]))
            surfaceDir[ix, iy] .+= 1
            surfaceDir[ix, iy] ./= 2.1
        end
    end
end
#

pictureNoise = [perlin((1 .*s .+ 11)...) for s in surfaceDir]
# p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(2*res*20,2*res*20), showaxis=false, widen=false, color=:sun)
save("back$(picNo).jld2", "yCenter", ycent, "xDir", xDir, "yDir", yDir, "picture", picture, "pictureNoise", pictureNoise)
p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(2*res,2*res), showaxis=false, widen=false, color=:turbid)

display(p)

save("candidate$(picNo).png",p)
# end

# save("picNoise$(counter).png", p)


## FRONT ################################

using AbstractAlgebra, Plots

function compPattern(m)

    s = factorial(m-1)

    indices = [vcat([1], v.d .+ 1) for v in collect(SymmetricGroup(m-1))]
    sort!(indices)

    M = zeros(Int, s,s)
    sInd1 = zeros(Int, m)
    sInd2 = zeros(Int, m)
    rev1 = zeros(Int, m)
    rev2 = zeros(Int, m)

    c = 0
    # for k = 1:length(indices)
    #     for i = 1:k
    #         ind1 = indices[i]
    #         for j = (i == k ? 1 : k) : k 
    #             ind2 = indices[j]

    for (i, ind1) in enumerate(indices)
        for (j, ind2) in enumerate(indices)
            if M[i,j] == 0
                @show (i,j, count(M .!= 0), count(M .== 0))
                c+= 1
                M[i,j] = c
                M[j,i] = c
                for sigma in SymmetricGroup(m)
                    sInd1 .= [sigma[k] for k in ind1]
                    sInd2 .= [sigma[k] for k in ind2]
                    sInd1 .= circshift(sInd1, -findfirst(x->x==1, sInd1)+1)
                    sInd2 .= circshift(sInd2, -findfirst(x->x==1, sInd2)+1)
                    x = findfirst(x->x==sInd1, indices)
                    y = findfirst(x->x==sInd2, indices)
                    @assert M[x,y] == c || M[x,y] == 0
                    M[x,y] = c
                    M[y,x] = c
                    
                    rev1 .= reverse(sInd1)
                    rev2 .= reverse(sInd2)
                    rev1 .= circshift(rev1, -findfirst(x->x==1, rev1)+1)
                    rev2 .= circshift(rev2, -findfirst(x->x==1, rev2)+1)
                    x = findfirst(x->x==rev1, indices)
                    y = findfirst(x->x==rev2, indices)
                    @assert M[x,y] == c || M[x,y] == 0
                    M[x,y] = c
                    M[y,x] = c
                end
            end
        end
        # end
    end
    return M
end


m = 8
M = compPattern(m)

reverse!(M, dims=1)

# p = heatmap(M , xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(factorial(m-1),factorial(m-1)), showaxis=false, widen=false)
# p = heatmap(M .* M , xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(factorial(m-1),factorial(m-1)), showaxis=false, widen=false)
# p = heatmap(M .^ 0.5 , xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(factorial(m-1),factorial(m-1)), showaxis=false, widen=false)
# p = heatmap(sin.(M / maximum(M)) , xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(factorial(m-1),factorial(m-1)), showaxis=false, widen=false)

scl = maximum(M)

# mulD = Dict(c=>count(M .== c) for c in unique(M))

mMod = M / scl

# mMod = log.(mMod)
# mMod = log.(mMod .^ 100)
# mMod = log.(log.(M) .+ 10)/ scl*100

mMod = maximum(mMod) .- mMod


res = factorial(m-1)

# mMod ./= 10

# MNoise = [octaveperlin(c+seed, c+seed, c+seed, 2,1) for c in mMod]
MNoise = mMod

p = heatmap(MNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(res,res), showaxis=false, widen=false, color=colorschemes[:YlOrBr][2:end-2])

display(p)

# save("front.jld2", "M", M)

##

using FileIO
# save("front.png",p)

##

using FileIO, Plots, LinearAlgebra
back = load("back8.jld2")
front = load("front.jld2")

picture = back["picture"]
M = front["M"]

##

surfaceDir = [[0.0,0.0,0.0] for i = 1:size(picture,1)-1,j=1:size(picture,1)-1]

res = Int((size(picture,1)-1)/2)

for (ix, x) in enumerate(-1:(1/res):1)
    for (iy, y) in enumerate(-1:(1/res):1)
        if ix < size(picture,1) && iy < size(picture,1)
            a = [x,y,picture[ix,iy]]
            b = [x+1/res,y,picture[ix+1,iy]]
            c = [x,y+1/res,picture[ix,iy+1]]
            surfaceDir[ix,iy] .= cross(b-a, c-a)
            surfaceDir[ix, iy] ./= sqrt(dot(surfaceDir[ix, iy],surfaceDir[ix, iy]))
            surfaceDir[ix, iy] .+= 1
            surfaceDir[ix, iy] ./= 2.1
        end
    end
end
##

using ColorSchemes

i = 1

cols = collect(colorschemes)
##

# :Oranges_7 :Oranges_4 :seaborn_icefire_gradient :tol_YlOrBr :YlOrRd_7 :turbid :flag_mk :linear_wyor_100_45_c55_n256 :flag_gr :flag_hn :tol_ylorbr :YlOrBr :ylOrBr_3 :Oranges_6 :OrRd

i+=1
#
@show cols[i][1]

# pictureNoise = [perlin((1 .*s .+ 11)...) for s in surfaceDir]

##
res = 1000
seed = rand()
pictureNoise = [octaveperlin((s .+ seed)...,2,2) for s in surfaceDir]

p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(2*res,2*res), showaxis=false, widen=false, color=colorschemes[:YlOrBr][2:end-2])


i += 1
save("back$i.png",p)

##

using ColorSchemes

seed = rand()*10
scl = maximum(M)

# mulD = Dict(c=>count(M .== c) for c in unique(M))

mMod = M / scl

# mMod = log.(mMod)
# mMod = log.(mMod .^ 100)
# mMod = log.(log.(M) .+ 10)/ scl*100

mMod = maximum(mMod) .- mMod


res = factorial(8-1)

# mMod ./= 10

# MNoise = [octaveperlin(c+seed, c+seed, c+seed, 2,1) for c in mMod]
MNoise = mMod
# 
p = heatmap(MNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(res,res), showaxis=false, widen=false, color=colorschemes[:YlOrBr][2:end-2])

save("front.png",p)

## Animated back cover 
include("perlin.jl")
##

#TODO: Need high precision center, opt

picNo = 8
using ColorSchemes, Plots, LinearAlgebra

res = 200
yRes = -1:(1/res):1
xRes = -(9/16):(1/res):(9/16)

data1 = load("pic8.jld2")
data2 = load("pic42.jld2")
data3 = load("pic39.jld2")
data4 = load("pic52.jld2")
xDir1 = data1["xDir"] * 1.2
yDir1 = data1["yDir"] * 1.2
xDir2 = data2["xDir"] * 1.2
yDir2 = data2["yDir"] * 1.2
xDir3 = data3["xDir"] * 1.2
yDir3 = data3["yDir"] * 1.2
xDir4 = data4["xDir"] * 1.2
yDir4 = data4["yDir"] * 1.2

ycent = data1["yCenter"]


# circle trough 1,2,3

A = vcat(vec(xDir1), vec(yDir1))
B = vcat(vec(xDir2), vec(yDir2))
C = vcat(vec(xDir3), vec(yDir3))

# u1 = B - A
# w1 = cross(C-A, u1)
# u = u1/norm(u1)
# w = w1/norm(w1)
# v = cross(w,y)

u1 = A-B
u = u1/norm(u1)
v1 = A-C 
v2 = v1 - dot(v1,u)*u
v = v2/norm(v2)

b = [dot(B-A, u),0]
c = [dot(C-A, u), dot(C-A,v)]

h = ((c[1]-b[1]/2)^2+c[2]^2-(b[1]/2)^2)/(2*c[2])

center = A + (b[1]/2)*u + h*v

radius = norm(center-A)

function circPoint(sigma)
    tmp = center + radius * sin(sigma)*u + radius * cos(sigma)*v
    X = reshape(tmp[1:9], 3,3)
    Y = reshape(tmp[10:end], 3,3)
    return X,Y
end

function compImage(yc, xD, yD, r)

    P = zeros(Float64, length(xRes), length(yRes))

    @showprogress for (i, x) in enumerate(xRes)
        for (j, y) in enumerate(yRes)
            P[i,j] = distanceTillFeasible(yc + x*xD + y*yD)
        end
    end

    # @time picture = [distanceTillFeasible(yc + x*xD + y*yD) for x=-1:(1/r):1, y=-1:(1/r):1]
    return P
end

frame = 1

@gif for sigma = 0:(2*pi/200):2*pi
    # if Int(round(10000*sigma/(2*pi))) % 100 == 0
    if isfile("./animFrames/res$(res)frame$(Int(round(10000*sigma/(2*pi)))).png")
        continue
    end
    
    @show sigma/(2*pi)

    # @time picture = compImage(ycent,
    #     ((1+sin(sigma))/2)*xDir1 + (1-((1+sin(sigma))/2))*xDir3 + ((1+cos(sigma))/2)*xDir2 + (1-((1+cos(sigma))/2)) * xDir4,
    #     ((1+sin(sigma))/2)*yDir1 + (1-((1+sin(sigma))/2))*yDir3 + ((1+cos(sigma))/2)*yDir2 + (1-((1+cos(sigma))/2)) * yDir4,
    #      res)
    # @time picture = compImage(ycent,
    #     sin(sigma)*xDir1 +cos(sigma)*xDir2,
    #     sin(sigma)*yDir1 + cos(sigma)*yDir2,
    #      res)
    @time picture = compImage(ycent,
        circPoint(sigma)...,
         res)
         

    # @time picture = [distanceTillFeasible(ycent + x*xDir + y*yDir) for x=-1:(1/res):1, y=-1:(1/res):1]

    using Plots

    pic = deepcopy(picture)
    pic[pic .<= -100] .= minimum(pic[pic .> -100])

    # heatmap(Float64.(pic))

    # distanceTillFeasible(0, ycent, 1, zero(ycent))

    #
    gr()
    # p = plot(Float64.(pic), st=:surface)


    # save("pic$(counter).png", p)

    surfaceDir = [[0.0,0.0,0.0] for i = 1:size(picture,1)-1,j=1:size(picture,2)-1]

    for (ix, x) in enumerate(xRes)
        for (iy, y) in enumerate(yRes)
            if ix < size(picture,1) && iy < size(picture,2)
                a = [x,y,picture[ix,iy]]
                b = [x+1/res,y,picture[ix+1,iy]]
                c = [x,y+1/res,picture[ix,iy+1]]
                surfaceDir[ix,iy] .= cross(b-a, c-a)
                surfaceDir[ix, iy] ./= sqrt(dot(surfaceDir[ix, iy],surfaceDir[ix, iy]))
                surfaceDir[ix, iy] .+= 1
                surfaceDir[ix, iy] ./= 2.1
            end
        end
    end
    #

    pictureNoise = [perlin((1 .*s .+ 11)...) for s in surfaceDir]
    # p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(2*res*20,2*res*20), showaxis=false, widen=false, color=:sun)
    # save("back$(picNo).jld2", "yCenter", ycent, "xDir", xDir, "yDir", yDir, "picture", picture, "pictureNoise", pictureNoise)
    p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(length(yRes),length(xRes)), showaxis=false, widen=false, color=colorschemes[:YlOrBr][2:end-2])

    display(p)
    save("./animFrames/res$(res)frame$(Int(round(10000*sigma/(2*pi)))).png", p)
    save("./animFrames/res$(res)frame$(Int(round(10000*sigma/(2*pi))))_data.jld2", "picture", picture, "pictureNoise", pictureNoise)
    # frame += 1

end

##

using FileIO, ImageShow

for (e,i) in enumerate(0:50:9950)
    # img = load("animFrames/res100frame$i.png")

    # cp("animFrames/res100frame$i.png", "animFrames/video/res100frame$e.png")
    img = load("animFrames/res$(res)frame$i.png")
    # save("animFrames/video/res100frame$e.png",img[1:end-1,1:end-1])
    save("animFrames/video/res$(res)frame$e.png",img[1:end-1,1:end-1])
    # run(`mogrify -blur 0.1 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 4 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 10 animFrames/video/res100frame$e.png`)
    # run(`mogrify -bilateral-blur 10 animFrames/video/res100frame$e.png`)
    run(`mogrify -kuwahara 1 animFrames/video/res$(res)frame$e.png`)
    # run(`mogrify -kuwahara 1 animFrames/video/res100frame$e.png`)
    # run(`mogrify -kuwahara 1 animFrames/video/res100frame$e.png`)
    # run(`mogrify -kuwahara 1 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 2 animFrames/video/res100frame$e.png`)
    run(`mogrify -blur 1x2 animFrames/video/res$(res)frame$e.png`)
    # run(`mogrify -filter Gaussian -resize 400% animFrames/video/res100frame$e.png`)
    run(`mogrify -filter Mitchell -adaptive-resize 400% animFrames/video/res$(res)frame$e.png`)
    # run(`mogrify -unsharp 3 animFrames/video/res100frame$e.png`)
    # run(`mogrify -unsharp 10 animFrames/video/res100frame$e.png`)
    # run(`mogrify -unsharp 10 animFrames/video/res100frame$e.png`)
    # run(`mogrify -kuwahara 1 animFrames/video/res100frame$e.png`)
    # run(`mogrify -blur 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -blur 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -blur 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -bilateral-blur 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 5 animFrames/video/res100frame$e.png`)
    # run(`mogrify -adaptive-sharpen 5 animFrames/video/res100frame$e.png`)
end

run(`ffmpeg -framerate 10 -i ./animFrames/video/res$(res)frame%00d.png -y -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p output.mp4`)

