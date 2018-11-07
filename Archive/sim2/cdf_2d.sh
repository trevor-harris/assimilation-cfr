#!/bin/bash

i=1
for f1 in {100,500,1000}; do
for f2 in {100,500}; do
for d in {10,20,30}; do
for l in {1,5,10}; do 

qsub -v F1=$f1,F2=$f2,D=$d,L=$l,I=$((i++)) cdf_2d.pbs

done
done
done
done


