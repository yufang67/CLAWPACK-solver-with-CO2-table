#!/bin/bash
#SBATCH --job-name=noz1
#SBATCH --mail-user=owngoals@163.com
#SBATCH --mail-type=ALL
#SBATCH --time=72:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue


srun xclaw > stdout.$SLURM_JOB_ID
