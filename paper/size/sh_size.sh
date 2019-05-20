#!/bin/bash

# simulate with and without phase for each of the 8 models
i=1

for n1 in 50 100 200 300 400 500; do
for n2 in 50 100 200 300 400 500; do
for rng in 0.25 0.75 1.25 1.75 2.25; do
for nu in 0.25 0.75 1.25 1.75 2.25; do

qsub -v N1=$n1,N2=$n2,RNG=$rng,NU=$nu,I=$((i++)) pbs_size.pbs

done
done
done
done
