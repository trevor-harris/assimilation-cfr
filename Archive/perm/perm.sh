#!/bin/bash

i=1
for i in {1..20}; do
qsub -v I=$i perm.pbs
done
