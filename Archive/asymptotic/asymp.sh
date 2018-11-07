#!/bin/bash

i=1
for i in {1..50}; do
qsub -v I=$i asymp.pbs
done
