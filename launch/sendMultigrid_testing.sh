#!/bin/bash
#SBATCH -N 1 
#SBATCH --ntasks-per-node=16
#SBATCH -t 48:00:00
#SBATCH -p xula4

#SBATCH -o output_mgTest
#SBATCH -e error_mgTest
#SBATCH --job-name=mgTest


SLURM_NPROCS=`expr $SLURM_JOB_NUM_NODES \* $SLURM_NTASKS_PER_NODE`
N_GROUPS=`expr $SLURM_NPROCS / 8` #Number of parallel cases to run (8 per process)

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

#Parallel julia execution. Script arguments: 1) Ranks per multigrid level, 2) Number of parallel cases, 3) Id of the parallel case



pids=()

for id in $(seq 1 "$N_GROUPS"); do
    echo "Launching multigrid testing for group $id"

    srun --exclusive --ntasks=8 --cpu-bind=cores \
        --output="mg_group${id}_%j.out" \
        --error="mg_group${id}_%j.err" \
        julia -O3 --check-bounds=no \
        ../scripts/multigrid_testing.jl "2,2,2" "$N_GROUPS" "$id" &

    pids+=("$!")
done

STATUS=0
for pid in "${pids[@]}"; do
    wait "$pid" || STATUS=1
done

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