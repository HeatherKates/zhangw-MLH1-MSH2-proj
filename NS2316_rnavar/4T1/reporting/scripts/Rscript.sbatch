#!/bin/bash
#SBATCH --job-name=rscript
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=12:00:00
#SBATCH --output=logs/rscript_%j.log
#SBATCH --error=logs/rscript_%j.error
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett

module purge; module load R/4.5
Rscript 01_annotate_and_extract_variants.R
