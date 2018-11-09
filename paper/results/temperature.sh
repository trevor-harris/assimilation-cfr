#!/bin/bash

# how many jobs to split out into
n=10

for t in $(seq 1 $n); do
qsub -v ERA=$t,N=$n temperature.pbs
done
