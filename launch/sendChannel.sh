#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=32
#SBATCH -t 48:00:00
#SBATCH -p xula4

#SBATCH -o output_channel
#SBATCH -e error_channel
#SBATCH --job-name=Channel


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`

srun hostname -s > hosts.$SLURM_JOB_ID
echo "================================================================"
hostname
echo "Using: ${SLURM_NPROCS} procs in ${SLURM_JOB_NUM_NODES} nodes"
echo "Running ${N_GROUPS} parallel cases"
echo "================================================================"
echo ""


SECONDS=0
#Load the enviroment
source ../env.sh

#Multigrid case
#    srun julia -O3 --check-bounds=no  ../scripts/Run_mg_channel.jl "4,4,2" 

# Direct solver with MUMPS case
     srun julia -O3 --check-bounds=no  ../scripts/Run_mumps_channel.jl "4,4,2"

STATUS=0

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
