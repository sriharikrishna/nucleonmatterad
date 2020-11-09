#!/bin/bash

#SBATCH --job-name=bfgs_pnm_2
#SBATCH --account=qocjax
#SBATCH --partition=bdwall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --output=bfgs_pnm_2.out
#SBATCH --error=bfgs_pnm_2.error
#SBATCH --mail-user=snarayan@.anl.gov
#SBATCH --time=20:00:00

#
module load intel-parallel-studio/cluster.2017.4-wyg4gfu intel/17.0.4-74uvhji
module load jdk/8u141-b15-mopj6qr

# Run My Program
echo "srun ./script_loop_pnm_fullx.sh ${START} ${STOP}"
srun ./script_loop_pnm_fullx.sh ${START} ${STOP}
