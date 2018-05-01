#!/bin/bash

for BATCHS in 1 2 3 4 5 6 7 8 9 10; do
	qsub -v BATCH=$BATCHS,SIMS=100 depth_bonf.pbs
done
