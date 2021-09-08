#!/bin/bash

#SBATCH --array=[1-10000]
#SBATCH -p scavenger
#SBATCH --mem=800M

/opt/apps/matlabR2021a/bin/matlab -nojvm -nodisplay -singleCompThread -r "rank=$SLURM_ARRAY_TASK_ID;BranchingColonyMultispecies_cluster;quit"

