#!/bin/bash

# reproducibility
i=1
s=042696

# parameters for F
n=300
pts=20

mu1=0
sd1=1
r1=30

# For job naming
mp="mean"
sp="sd"
rp="corr"

# Mean Power
for mu2 in $(seq -1 0.1 1); do
qsub -v NAME=$mp,N=$n,PTS=$pts,MU1=$mu1,MU2=$mu2,SD1=$sd1,SD2=$sd1,R1=$r1,R2=$r1,S=$s,I=$((i++)) power.pbs 
# echo $mp
done

# Standard Deviation Power
for sd2 in $(seq 0.1 0.05 2); do
qsub -v NAME=$sp,N=$n,PTS=$pts,MU1=$mu1,MU2=$mu1,SD1=$sd1,SD2=$sd2,R1=$r1,R2=$r1,S=$s,I=$((i++)) power.pbs
# echo $sd2
done


# Correlation Power
for r2 in $(seq 5 5 80); do
qsub -v NAME=$rp,N=$n,PTS=$pts,MU1=$mu1,MU2=$mu1,SD1=$sd1,SD2=$sd1,R1=$r1,R2=$r2,S=$s,I=$((i++)) power.pbs
# echo $r2
done
