#!/bin/bash

for BATCHS in {1..50}; do
	qsub -v BATCH=$BATCHS,SIMS=20 realistic.pbs
done
