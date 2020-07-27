#!/bin/bash
for ((i=${1};i<=${2};i++))
do
   for j in 0.1 0.2
   do
     #echo "RUNNING ./script_bfgs_pnm_fullx.sh $i $j 1.0 2.0 0.0 1.0 1.0 2.0 3.0 4.5"
     #./script_bfgs_pnm_fullx.sh $i $j 1.0 2.0 0.0 1.0 1.0 2.0 3.0 4.5
     echo "RUNNING ./script_bfgs_pnm_fullx.sh $i $j 1.0 2.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0 2.0 3.0 4.5 1.0 2.0"
     ./script_bfgs_pnm_fullx.sh $i $j 1.0 2.0 0.0 1.0 0.0 1.0 0.0 1.0 1.0 2.0 3.0 4.5 1.0 2.0
   done
done

