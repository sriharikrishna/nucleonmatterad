#!/bin/bash
for ((i=${1};i<=${2};i++))
do
   for j in 0.1 0.2
   do
     echo "RUNNING ./script_bfgs_snm_fullx.sh $i $j 3.0 4.8 -0.5 1.9 -0.5 1.9 -0.5 1.9 1.0 1.0 1.0 1.0 0.0 0.0"
     ./script_bfgs_snm_fullx.sh $i $j 3.0 4.8 -0.5 1.9 -0.5 1.9 -0.5 1.9 1.0 1.0 1.0 1.0 0.0 0.0
   done
done

