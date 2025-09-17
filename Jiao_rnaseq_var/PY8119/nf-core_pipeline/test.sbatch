#!/bin/bash
#SBATCH --job-name=checkqos
#SBATCH --output=checkqos.out
#SBATCH --ntasks=1
#SBATCH --time=00:01:00

scontrol show job \$SLURM_JOB_ID
