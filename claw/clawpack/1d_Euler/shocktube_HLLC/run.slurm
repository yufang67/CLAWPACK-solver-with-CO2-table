#!/bin/bash
#SBATCH --job-name=clawpack_test
#SBATCH --mail-user=yu.fang@uclouvain.be
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --no-requeue


srun xclaw > stdout.$SLURM_JOB_ID
