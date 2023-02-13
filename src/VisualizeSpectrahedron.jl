using MosekTools, JuMP, LinearAlgebra

## random SDP: max <C,X> s.t. A(X) <= b, X PSD

n = 3
m = 30

C = Symmetric(rand(n,n))
Ab = [
    (Symmetric(rand(n,n)), rand()) for i = 1:m
    ]

########### 1) Find optimum xOpt, und center xCenter (interior point without objective)

m = Model(Mosek.Optimizer)
X = @variable(m, X[1:n, 1:n], PSD)
for (A,b) in Ab
    @constraint(m, dot(X, A) <= b)
end

optimize!(m)
xCenter = value(X)

@objective(m, Max, dot(X, C))

optimize!(m) 
xOpt = value(X)

########### 2) Find a random hyperplane orthogonal to xOpt - xCenter

centerOpt = xOpt - xCenter

xDir = Symmetric(2 .*rand(n, n) .- 1)
xDir -= centerOpt*dot(xDir, centerOpt)/dot(centerOpt, centerOpt)

yDir = Symmetric(2 .*rand(n, n) .- 1)
yDir -= centerOpt*dot(yDir, centerOpt)/dot(centerOpt, centerOpt)
yDir -= xDir * dot(yDir, xDir) / dot(xDir, xDir)

# normalize and scale appropriately for the SDP (could automate, but am lazy)
xDir = 1/3*norm(centerOpt) * xDir ./ dot(xDir, xDir)
yDir = 1/3*norm(centerOpt) * yDir ./ dot(yDir, yDir)

########### 4) Solve univariate SDP for every pixel in the plane given by xDir, yDir
# Here we substitute X = xCenter + shift + t*centerOpt,
# where shift lies in the plane [-1,1]*xDir + [-1,1]*yDir + xCenter.
# t is the only variable, and gives the "height" of the pixel 

# Prepare model in advance, we then "fix" appropriate variables

m2 = Model(Mosek.Optimizer)
set_silent(m2)

t = @variable(m2,t)

startY = @variable(m2, startY[1:n, 1:n], Symmetric) # will be fixed appropriately

X = startY + t * centerOpt
@constraint(m2, X in PSDCone())

for (A,b) in Ab
    @constraint(m2, dot(X, A) <= b)
end

@objective(m2, Max, t)



res = 100

heights = map(Iterators.product(-1:(1/(res/2-0.5)):1+(1/(res/2-0.5)), -1:(1/(res/2-0.5)):1+(1/(res/2-0.5)))) do (x,y)
    fix.(startY, Symmetric(Float64.(xCenter + x*xDir + y*yDir)), force=true)

    print("$y               \r")

    optimize!(m2)

    if termination_status(m2) in [JuMP.MathOptInterface.INFEASIBLE]#,JuMP.MathOptInterface.SLOW_PROGRESS]
        return -100
    end
 
    return value(t)
end

########### 4) Drawing the result
using Plots, ColorSchemes
include("perlin.jl")

# Easy, but not very readable: Use height directly
display(heatmap(heights))

sleep(1)

# Better: Compute surface vectors, apply perlin noise    
surfaceDir = map(Iterators.product(enumerate(-1:(1/(res/2-0.5)):1), enumerate(-1:(1/(res/2-0.5)):1))) do ((ix,x),(iy,y))
    fix.(startY, Symmetric(Float64.(xCenter + x*xDir + y*yDir)), force=true)

    a = [x,y,heights[ix,iy]]
    b = [x+1/res,y,heights[ix+1,iy]]
    c = [x,y+1/res,heights[ix,iy+1]]
    v = cross(b-a, c-a)
    v ./= sqrt(dot(v,v))
    v .+= 1 # arbitrary shift -- shifts noise
    v ./= 2.1 # arbitrary scale -- changes noise "speed"
    v
end

pictureNoise = [perlin(s...) for s in surfaceDir]
p = heatmap(pictureNoise, xaxis=:none, legend=:none, border=:none, margin=Plots.Measures.pt * -100, size=(4*res,4*res), showaxis=false, widen=false,color=colorschemes[:YlOrBr][2:end-2])

########### 5) Finding a nice camera angle, noise shifts, noise scale, color scheme... I generated many small pictures and picked what looked interesting. I applied some post-processing (anti-aliasing, etc), as edges are necessarily noisy sadly


