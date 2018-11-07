#!/bin/bash

i=1
s=1023
for n in {20,100,300,500}; do
for l in {30,40,50}; do

qsub -v N=$n,D=3,L=$l,S=$s,I=$((i++)) size_rev.pbs

# echo $n $d $l $s $((i++))

done
done


