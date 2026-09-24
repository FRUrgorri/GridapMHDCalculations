#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=8
#SBATCH -t 48:00:00
#SBATCH -p xula4

#SBATCH -o output_job
#SBATCH -e error_job
#SBATCH --job-name=name
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

srun julia --project=.. -O3 --check-bounds=no ../test/mpi/runtest_mpi.jl

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

