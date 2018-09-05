#!/bin/bash

for SHIFTS in 0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1; do
	qsub -v SHIFT=$SHIFTS,SIMS=100 mu.pbs
done
