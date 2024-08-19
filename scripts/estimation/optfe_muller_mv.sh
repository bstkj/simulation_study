#!/bin/bash

k=$1 # max cluster size
# split up jobs into batches, so that each physical core is only assigned 1 job at a time
# (for reduced runtimes)
batch=$2
echo "STARTED: k=$k, batch=$batch"
# path to Julia executable
juliapath=/workspace/bteo/julia-1.10.5/bin/julia
for j in {1..20}
do
    i=$((batch * 20 + j)) # dataset index
    echo "j=$j, i=$i"
    # taskset -c $((j + 3)): assigns job to logical core $((j + 3))
    # 2>&1: redirect stderr to stdout
    # | tee "... .log": write stdout to file, but still display it
    # &: runs command in background
    taskset -c $((j + 3)) $juliapath $(pwd)/scripts/estimation/optfe_muller_mv.jl $i $k 2>&1 | tee "$(pwd)/results/estimation/muller_mv/rep$i.log" &
done
wait # wait for background jobs started in loop to finish
echo "COMPLETED: k=$k, batch=$batch"