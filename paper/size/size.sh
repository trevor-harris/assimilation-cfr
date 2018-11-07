#!/bin/bash

i=1
s=0426
d=1

for n in {50,100,300}; do
for l in {20,30,40,50}; do

qsub -v N=$n,D=$d,L=$l,S=$s,I=$((i++)) size.pbs

# echo $n $d $l $s $((i++))

done
done


