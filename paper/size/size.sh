#!/bin/bash

i=1
for s in {66220,84710,93477,69707,95547,4581,53132,64050,70638,51740}; do
for n in {100,250,500}; do
for l in {5,10,20}; do
for d in {5,5,5}; do 

qsub -v N=$n,D=$d,L=$l,S=$s,I=$((i++)) size.pbs

done
done
done
done


