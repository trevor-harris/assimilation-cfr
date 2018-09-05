#!/bin/bash

for BATCHES in 1..100; do
	qsub -v BATCH=$BATCHES actual.pbs
done
