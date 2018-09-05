#!/bin/bash

for SCALES in 0.01 0.05 0.1 0.3 0.5 0.7 1.3 1.5 1.7 2 3 4 5; do
	qsub -v SCALE=$SCALES,SIMS=100 cov.pbs
done
