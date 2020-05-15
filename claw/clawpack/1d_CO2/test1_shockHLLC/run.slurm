#!/bin/bash
#SBATCH --job-name=clawpack_test
#SBATCH --mail-user=owngoals@163.com
#SBATCH --mail-type=ALL
#SBATCH --time=00:20:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue


srun xclaw > stdout.$SLURM_JOB_ID
