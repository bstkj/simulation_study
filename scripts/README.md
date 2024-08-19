This folder contains scripts used for simulating data, fitting models, and visualizing
results.

1. `simtraits.jl`: simulate trait data (`data/traits`) on networks (`data/networks`)
2. `pick_k/`: profile accuracy and runtime vs maximum cluster size for computing the
factored energy
    - univariate trait: `pickk_sikora_uv.jl`, `pickk_lipson_uv.jl`, `pickk_muller_uv.jl`
    - multivariate trait: `pickk_sikora_mv.jl`, `pickk_lipson_mv.jl`, `pickk_muller_mv.jl`
    - run analyses for all datasets in batches: `pickk.sh`
    - aggregate results across datasets: `aggregate.R`
    - make figures (`figures/pickk*.png`): `pickk.R`
3. `estimation/`: compare maximum factored energy (MFE) estimation with maximum likelihood
(ML) estimation, for the Müller network in the multivariate case
    - MFE estimation: `optfe_muller_mv.jl`
    - run MFE analyses for all datasets in batches: `optfe_muller_mv.sh`
    - ML estimation: `optll_muller_mv.jl`
    - aggregate results across datasets: `aggregate.jl`
    - make figures (`figures/opt_muller_mv*.png`): `rms_deviations.R`, `residuals.R`

*`pickk.sh` and `optfe_muller_mv.sh`, which are for batch analyses, may need to be
customized to the user's computing resources