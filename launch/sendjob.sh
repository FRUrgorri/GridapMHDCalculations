#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH -p xula4

#SBATCH -o output_channel
#SBATCH -e error_channel
###SBATCH --mail-user=fernando.roca@ciemat.es
#SBATCH --job-name=channel
###SBATCH --mem=0

SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs in ${SLURM_JOB_NUM_NODES} nodes"
echo "================================================================"
echo ""


SECONDS=0
#Load the enviroment
source ../env.sh

#Parallel julia execution

PASS_FILE="Transfer.jl"

echo "_np = (4,4,1)" >> $PASS_FILE
echo "@assert isequal(_np[1]*_np[2]*_np[3],$SLURM_NPROCS)" >> $PASS_FILE

srun julia --project=.. -O3 --check-bounds=no ./scripts/channel_Ha10Re1_H1Hdiv_direct.jl

rm $PASS_FILE

duration=$SECONDS
rm -f hosts.$SLURM_JOB_ID
rm -f input_params.jl
rm -f mesh_params

STATUS=$?
echo "================================================================"
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
echo "================================================================"
echo ""
echo "STATUS = $STATUS"
echo ""

