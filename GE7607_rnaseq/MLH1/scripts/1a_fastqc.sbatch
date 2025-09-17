#!/bin/bash
#SBATCH --job-name=fastqc.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8GB
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH -t 02:00:00
#SBATCH --output=logs/fastqc.%j.out
#SBATCH --error=logs/fastqc.%j.err
#SBATCH --array=1-6

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file performs FastQC quality control for MLH1 data samples
###################################################################

module load fastqc/0.12.1

# Define the base directory for input data
DATA_DIR="data"

# Define the output directory for FastQC results
OUT_DIR="results/1_fastqc"
mkdir -p ${OUT_DIR}

# Create an array of R1 files in the data directory
R1_FILES=(${DATA_DIR}/*_R1_001.fastq.gz)

# Get the file for the current task based on SLURM_ARRAY_TASK_ID
Sample1=${R1_FILES[$((SLURM_ARRAY_TASK_ID-1))]}
Sample2=$(echo ${Sample1} | sed 's/_R1_/_R2_/')

# Extract sample name from the input file path
SAMPLE_NAME=$(basename ${Sample1} | sed 's/_R1_001.fastq.gz//')

# Run FastQC on R1 and R2 files
fastqc ${Sample1} ${Sample2} -o ${OUT_DIR} -t 2

echo -e "\nInfo: FastQC job completed for ${SAMPLE_NAME} on $(date) on $(hostname).\n"

