#!/bin/bash

k=$1 # max cluster size
# split up jobs into batches (for benchmarking), so that each physical core is only assigned
# 1 job at a time
batch=$2
echo "STARTED: k=$k, batch=$batch"
# path to Julia executable
juliapath=/workspace/bteo/julia-1.10.5/bin/julia
for j in {1..20}
do
    i=$((batch * 20 + j)) # dataset index
    echo "j=$j, i=$i"
    # taskset -c $j: assigns job to logical core $j (useful for benchmarking)
    # 2>&1: redirect stderr to stdout
    # | tee "... .log": write stdout to file, but still display it
    # &: runs command in background
    taskset -c $j $juliapath $(pwd)/scripts/pick_k/pickk_sikora_uv.jl $i $k 2>&1 | tee "$(pwd)/results/pick_k/sikora_uv/k$k/rep$i.log" &
    # change pickk_sikora_uv.jl and sikora_uv above as needed for different networks, trait
    # dimension, e.g., pickk_lipson_mv.jl and lipson_mv
done
wait # wait for background jobs started in loop to finish
echo "COMPLETED: k=$k, batch=$batch"