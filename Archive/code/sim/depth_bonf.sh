#!/bin/bash

for RUNS in 1 2 3; do
for BATCHS in {1..25}; do
	qsub -v BATCH=$BATCHS,SIMS=100,RUN=$RUNS depth_bonf.pbs
done
done
