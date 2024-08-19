using CSV
using DataFrames
using PhyloNetworks
const PN = PhyloNetworks
using Random

net_sikora = readTopology("data/networks/sikora_2019.phy");
net_lipson = readTopology("data/networks/lipson_2020b.phy");
net_muller = readTopology("data/networks/muller_2022.phy");

B = 100 # no. of datasets per simulation setting

# Simulate univariate data
for (net, id) in [(net_sikora, "sikora"), (net_lipson, "lipson"), (net_muller, "muller")]
    Random.seed!(17920921)
    traits = zeros(net.numTaxa, B)
    for i in 1:B
        traits[:,i] = simulate(net, ParamsBM(0,1))[:Tips]
    end
    df = hcat(DataFrame(:tipNames => tipLabels(net)), # rownames
              DataFrame(traits, [Symbol("d$i") for i in 1:B]))
    CSV.write("data/traits/"*id*"_univariatedata", df) # save dataset
end

# Simulate multivariate data
rootmean = zeros(4)
rate = [0.8 -0.71 -0.8 0.49;
        -0.71 0.8 0.81 -0.41;
        -0.8 0.81 1.1 -0.4;
        0.49 -0.41 -0.4 0.5]
for (net, id) in [(net_sikora, "sikora"), (net_lipson, "lipson"), (net_muller, "muller")]
    Random.seed!(353)
    tipNames = tipLabels(net)
    for i in 1:B
        sim = simulate(net, ParamsMultiBM(rootmean, rate))
        df = hcat(DataFrame(:tipNames => tipNames),
                  DataFrame(transpose(sim[:Tips]), :auto))
        CSV.write("data/traits/"*id*"_multivariatedata/rep_"*string(i), df)
    end
end