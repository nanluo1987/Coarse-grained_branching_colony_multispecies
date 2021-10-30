#!/bin/bash

#SBATCH --array=[1-14000]
#SBATCH -p scavenger
#SBATCH --mem=2G

/opt/apps/rhel8/matlabR2020b/bin/matlab -nojvm -nodisplay -singleCompThread -r "rank=$SLURM_ARRAY_TASK_ID;cluster_localScreening;quit"

