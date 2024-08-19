using CSV, DataFrames
using LinearAlgebra
using PhyloGaussianBeliefProp
using PhyloNetworks
const PGBP = PhyloGaussianBeliefProp
using Statistics
using Tables
using Optim

i = parse(Int, ARGS[1]) # max no. of iterations before FE is computed
k = parse(Int, ARGS[2]) # max cluster size

net = readTopology("real_networks/muller_2022.phy")
params = CSV.read("results/estimation/muller_mv_ml.csv", DataFrame)
niter = 50
cg = PGBP.clustergraph!(net, PGBP.JoinGraphStructuring(k))
df = CSV.read("data/muller_multivariatedata/rep_$i", DataFrame)
tbl_var = columntable(select(df, Not(:tipNames)))
μhat, Rhat = eval(Meta.parse(params.Θhat[i]))
μhat = vec(μhat)
m = PGBP.MvFullBrownianMotion(Rhat, μhat)
b, (n2c, n2fam, n2fix, n2d, c2n) = PGBP.allocatebeliefs(tbl_var, df.tipNames,
                                                       net.nodes_changed, cg, m)
cgb = PGBP.ClusterGraphBelief(b, n2c, n2fam, n2fix, c2n)

# max factored energy estimate using cluster graph
fitted, fenergy, opt = PGBP.calibrate_optimize_clustergraph!(cgb, cg, net.nodes_changed,
                        tbl_var, df.tipNames,
                        PGBP.MvFullBrownianMotion, (Symmetric(Rhat), μhat), niter,
                        PGBP.regularizebeliefs_onschedule!,
                        Optim.Options(f_tol=0.0001, iterations=50, show_trace=true, show_every=1))
CSV.write("results/estimation/muller_mv/rep_$i.csv",
          DataFrame(μ=[fitted.μ], R=[fitted.R], fenergy=[fenergy],
                    stopped_by=[Tuple(opt.stopped_by)], fg_calls=[(opt.f_calls, opt.g_calls)],
                    time_run=[opt.time_run])
         )