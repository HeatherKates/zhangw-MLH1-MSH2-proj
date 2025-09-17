#!/bin/bash
#SBATCH --job-name=genrich_array
#SBATCH --output=../logs/genrich_%A_%a.out
#SBATCH --error=../logs/genrich_%A_%a.err
#SBATCH --array=0-5  # For 6 BAM files (index starts at 0)
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  # Adjust based on your resource needs
#SBATCH --mem=250GB  # Adjust memory as needed
#SBATCH --time=02:00:00  # Adjust based on expected runtime
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett-b

# Load necessary modules
module load genrich

# Get all sorted BAM files in the directory
bam_path="results/4a_samtools"
bam_files=($(ls $bam_path/*.sorted.bam))

# Get the BAM file for the current SLURM_ARRAY_TASK_ID
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID]}

# Extract the sample name from the BAM file name
sample=$(basename $bam_file .sorted.bam)

# Run Genrich on the current BAM file
Genrich -t "$bam_file" -o "results/5_genrich/${sample}.narrowPeak" -j -y -r -v -e MT -d 100 -q 0.05 

