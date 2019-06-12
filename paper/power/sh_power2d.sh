#!/bin/bash

# reproducibility
i=1
s=042696

# parameters for F
n1=100
n2=50
pts=20

mu1=0
sd1=1
l=0.4
nu=1.0

# For job naming
mp="mean"
sp="sd"
rp="corr"

# Mean Power
for mu2 in $(seq -1 0.1 1); do
qsub -v NAME=$mp,N1=$n1,N2=$n2,PTS=$pts,MU1=$mu1,MU2=$mu2,SD1=$sd1,SD2=$sd1,R1=$l,R2=$l,NU1=$nu,NU2=$nu,S=$s,I=$((i++)) pbs_power_const.pbs
# echo $mp
done

# Standard Deviation Power
for sd2 in $(seq 0.1 0.05 2); do
qsub -v NAME=$sp,N1=$n1,N2=$n2,PTS=$pts,MU1=$mu1,MU2=$mu1,SD1=$sd1,SD2=$sd2,R1=$l,R2=$l,NU1=$nu,NU2=$nu,S=$s,I=$((i++)) pbs_power_const.pbs
# echo $sd2
done

# Correlation Power
for r2 in $(seq 0.05 0.05 1); do
qsub -v NAME=$rp,N1=$n1,N2=$n2,PTS=$pts,MU1=$mu1,MU2=$mu1,SD1=$sd1,SD2=$sd1,R1=$l,R2=$r2,NU1=$nu,NU2=$nu,S=$s,I=$((i++)) pbs_power_const.pbs
# echo $r2
done

# Smoothness Power
for nu2 in $(seq 0.1 0.1 2); do
qsub -v NAME=$rp,N1=$n1,N2=$n2,PTS=$pts,MU1=$mu1,MU2=$mu2,SD1=$sd1,SD2=$sd1,R1=$l,R2=$l,NU1=$nu,NU2=$nu2,S=$s,I=$((i++)) pbs_power_const.pbs

