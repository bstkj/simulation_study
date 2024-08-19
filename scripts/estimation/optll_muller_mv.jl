using CSV
using DataFrames
using LinearAlgebra
using PhyloNetworks
using Statistics

net = readTopology("data/networks/muller_2022.phy")
V = sharedPathMatrix(net)[:Tips]
h = median(diag(V))
LLopt = []
Θopt = []
Θhat = []
μtrue = zeros(4)
Rtrue = [0.8 -0.71 -0.8 0.49;
        -0.71 0.8 0.81 -0.41;
        -0.8 0.81 1.1 -0.4;
        0.49 -0.41 -0.4 0.5]
for i in 1:100
    df = CSV.read("data/traits/muller_multivariatedata/rep_$i", DataFrame)
    X = Matrix(select(df, Not(:tipNames)))
    μhat = mean(X, dims=1)[:]
    Rhat = cov(X) / h
    push!(Θhat, (μ=μhat, R=Rhat))
    μopt = transpose(ones(1,40) * (V \ X)) / dot(ones(40), V \ ones(40))
    residopt = X - ones(40,1) * transpose(μopt)
    Ropt = (transpose(residopt) * (V \ residopt)) / 40
    push!(Θopt, (μ=μopt, R=Ropt))
    residopt = transpose(residopt)[:]
    llopt = -0.5*transpose(residopt)*inv(kron(V, Ropt))*residopt - 0.5*logdet(2π*kron(V, Ropt))
    push!(LLopt, llopt)
end
optll_muller_mv = DataFrame(LLopt=LLopt, Θopt=Θopt, Θhat=Θhat)
CSV.write("results/estimation/muller_mv_ml.csv", optll_muller_mv)

using RCall
R"""
saveRDS($optll_muller_mv, file="results/estimation/muller_mv_ml.rds")
"""