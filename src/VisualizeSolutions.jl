# Compute beta_m for up to m = mMax and m odd, and visualize the eigenvector of the rank one solutions.

include("SolveSDP_Beta_Iterative.jl")

using GenericLinearAlgebra, JLD2

Ys = Dict()

mMax = 11

for m = 5:2:mMax
    optVal, Y, cs = solveBetaIterative(m);
    Ys[m] = Y .* factorial(m-1)
end

for m = 5:2:mMax
    @show m
    display(eigen(Hermitian(Ys[m])).values)
end

##

data = Dict()
evs = Dict()

for m = 5:2:mMax
    ev = eigen(Hermitian(Ys[m]))
    data[m] = (sqrt(ev.values[end]) * ev.vectors[:, end] ./ sign(ev.vectors[1,end]))
    evs[m] = ev.values[end]
end

##

using Plots
gr()

p = plot(data[5], legend=:none, color=:blue, ylims = (0,2.5), grid=:none)
p = scatter!(p,data[5], legend=:none, color=:blue)

for m = 7:2:mMax
    ev = eigen(Hermitian(Ys[m]))
    plot!(p,data[m], color=:blue)
    scatter!(p,data[m], color=:blue)
end

xticks!([1:1:floor(mMax/2);], ["M_$(Int(i+2))" for i = 1:floor(mMax/2)])

annotate!([(1.5, 0.55, ("\$\\beta_5\$", 10, :blue, :left))])
annotate!([(2.5, 0.7, ("\$\\beta_7\$", 10, :blue, :left))])
annotate!([(3.5, 0.87, ("\$\\beta_9\$", 10, :blue, :left))])
annotate!([(4.5, 0.92, ("\$\\beta_{11}\$", 10, :blue, :left))])
# annotate!([(5.5, 1, ("\$\\beta_{13}\$", 10, :blue, :left))])

p
