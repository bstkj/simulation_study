This folder contains the results of fitting models on simulated data.

1. `pick_k/`: accuracy and runtime vs maximum cluster size for computing the factored energy
    - results for each dataset and maximum cluster size: 
    `sikora_uv/`, `sikora_mv/`, `lipson_uv/`, `lipson_mv/`, `muller_uv/`, `muller_mv/`
    - aggregated results across datasets: `sikora_uv.csv`, `sikora_mv.csv`, `lipson_uv.csv`,
    `lipson_mv.csv`, `muller_uv.csv`, `muller_mv.csv`
    - aggregated results across datasets, averaged for each maximum cluster size:
    `sikora_uv_mean.csv`, `sikora_mv_mean.csv`, `lipson_uv_mean.csv`,
    `lipson_mv_mean.csv`, `muller_uv_mean.csv`, `muller_mv_mean.csv`
2. `estimation/`: maximum factored energy (MFE) estimation and maximum likelihood (ML)
estimation, for the Müller network in the multivariate case
    - MFE and parameter estimates: `muller_mv/`
    - MFE and parameter estimates, aggregated across datasets: `muller_mv_fe.csv`,
    `muller_mv_fe.rds`
    - maximum log-likelihood and parameter estimates, aggregated across datasets:
    `muller_mv_ml.csv`, `muller_mv_ml.rds`