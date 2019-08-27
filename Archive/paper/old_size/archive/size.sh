#!/bin/bash

i=1
s=072393
d=1

for j in {1..20}; do
for n in {50,100,200,300,400}; do
for l in {20,30,40,50}; do

qsub -v N=$n,D=$d,L=$l,S=$s,I=$((i++)),J=$j size.pbs

# echo $n $d $l $s $((i++))

done
done
done

