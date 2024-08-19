using CSV, DataFrames
using LinearAlgebra
using PhyloGaussianBeliefProp
using PhyloNetworks
const PGBP = PhyloGaussianBeliefProp
using Statistics
using Tables
using Optim

i = parse(Int, ARGS[1])
k = parse(Int, ARGS[2])

# factored energy
net = readTopology("real_networks/muller_2022.phy")
niter = 50 
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(k))
sch = PGBP.spanningtrees_clusterlist(cg, net.nodes_changed)
df = CSV.read("data/muller_multivariatedata/rep_$i", DataFrame)
tbl_var = columntable(select(df, Not(:tipNames)))
μtrue = zeros(4)
Σtrue = [0.8 -0.71 -0.8 0.49;
        -0.71 0.8 0.81 -0.41;
        -0.8 0.81 1.1 -0.4;
        0.49 -0.41 -0.4 0.5]
m = PGBP.MvFullBrownianMotion(Σtrue, μtrue)
b, (n2c, n2fam, n2fix, n2d, c2n) = PGBP.allocatebeliefs(tbl_var, df.tipNames,
                                                        net.nodes_changed, cg, m)
PGBP.assignfactors!(b, m, tbl_var, df.tipNames, net.nodes_changed, n2c, n2fam, n2fix)
cgb = PGBP.ClusterGraphBelief(b, n2c, n2fam, n2fix, c2n)
stats_fe_before = @timed begin
    PGBP.regularizebeliefs_onschedule!(cgb, cg)
    _, calibrated = PGBP.calibrate!(cgb, sch, niter; auto=true, info=true)
    fetrue = PGBP.factored_energy(cgb)[3]
end
PGBP.assignfactors!(b, m, tbl_var, df.tipNames, net.nodes_changed, n2c, n2fam, n2fix)
PGBP.init_factors_frombeliefs!(cgb.factor, cgb.belief)
PGBP.init_messagecalibrationflags_reset!(cgb, true)
stats_fe_after = @timed begin
    PGBP.regularizebeliefs_onschedule!(cgb, cg)
    _, calibrated = PGBP.calibrate!(cgb, sch, niter; auto=true, info=true)
    fetrue = PGBP.factored_energy(cgb)[3]
end
# log-likelihood
y = transpose(Matrix(select(df, Not(:tipNames))))[:] 
V = sharedPathMatrix(net)[:Tips]
V⊗Σtrue = kron(V, Σtrue)
μ𝟷ᵀ = kron(ones(length(net.leaf)), μtrue)
stats_ll_before = @timed begin
    lltrue = -0.5*transpose(y-μ𝟷ᵀ)*inv(V⊗Σtrue)*(y-μ𝟷ᵀ) - 0.5*logdet(2π*V⊗Σtrue)
end
stats_ll_after = @timed begin
    lltrue = -0.5*transpose(y-μ𝟷ᵀ)*inv(V⊗Σtrue)*(y-μ𝟷ᵀ) - 0.5*logdet(2π*V⊗Σtrue)
end
# residual
residtrue = lltrue - fetrue
# save results
CSV.write("data/pick_k/muller_mv/k$(k)/rep$(i).csv",
          DataFrame(fetrue=[fetrue], lltrue=[lltrue], residtre=[residtrue],
                    fetimeb=[stats_fe_before.time], fetimea=[stats_fe_after.time],
                    lltimeb=[stats_ll_before.time], lltimea=[stats_ll_after.time])
         )
