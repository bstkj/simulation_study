using CSV
using DataFrames

FEopt = []
Θopt = []
for i in 1:100
    df = CSV.read("results/estimation/muller_mv/rep$i.csv", DataFrame)
    push!(Θopt, (μ=eval(Meta.parse(df.μ[1])), R=eval(Meta.parse(df.R[1]))))
    push!(FEopt, df.fenergy[1])
end

optfe_muller_mv = DataFrame(FEopt=FEopt, Θopt=Θopt)
CSV.write("results/estimation/muller_mv_fe.csv", optfe_muller_mv)

using RCall
R"""
saveRDS($optfe_muller_mv, file="results/estimation/muller_mv_fe.rds")
"""