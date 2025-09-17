#!/bin/bash
#SBATCH --job-name=rserver
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100gb
#SBATCH --time=48:00:00
#SBATCH --output=rserver_%j.log
#SBATCH --error=rserver_%j.error
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett-b

module load R/4.2
rserver
