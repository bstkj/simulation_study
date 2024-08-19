# Adapting cluster graphs for inference of continuous trait evolution on phylogenetic networks

This repository contains code to reproduce simulations and figures in Chapter 4, "Adapting
cluster graphs for inference of continuous trait evolution on phylogenetic networks", of the
dissertation *Statistical and Computational Techniques for Continuous Trait Models on
Phylogenetic Networks*.

There are 4 top-level folders (`data`, `figures`, `results`, `scripts`), each containing
their own READMEs.

## How to run the scripts

All scripts should be run from the top-level of this repository (i.e. `simulation_study/`).

## Reproducing the same Julia environment

The `Project.toml` is provided so that the exact package environment used to perform the Julia-related analyses can be instantiated.

To reproduce the environment:
1. Open the Julia REPL within `simulation_study` (which contains `Project.toml`)
2. Enter the Pkg REPL and run `activate .`, then `instantiate`
```
(@v1.10) pkg> activate .
(simulation_study) pkg> instantiate
```

<!-- Simulation Study using [PhyloGaussianBeliefProp.jl](https://github.com/JuliaPhylo/PhyloGaussianBeliefProp.jl) to:
- quantify the accuracy vs scalability trade-off between exact and approximate inference
- understand how the adaptive elements of belief propagation can be designed to respond to the characteristics of the phylogeny -->