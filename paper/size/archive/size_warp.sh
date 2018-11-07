#!/bin/bash

i=1
for s in {6055,39158,2637,69903,11413,53004,31072,79767,16359,63731,61875,19749,36012,57035,1331,83542,47868,27346,87024,5358}; do
for n in {100,300,500}; do
for l in {50,90,130}; do
for d in {2,3,5}; do 

qsub -v N=$n,D=$d,L=$l,S=$s,I=$((i++)) size_warp.pbs

# echo $n $d $l $s $((i++))

done
done
done
done


