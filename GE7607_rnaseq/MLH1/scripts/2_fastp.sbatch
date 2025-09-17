#!/bin/bash
#SBATCH --job-name=fastp.%j
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8gb
#SBATCH --qos=timgarrett
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT
#SBATCH -t 01:00:00
#SBATCH --output=logs/fastp.%j.out
#SBATCH --error=logs/fastp.%j.err
#SBATCH --array=1-12

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform quality filtering and trimming with fastp
###################################################################

module load fastp/0.23.4

# Define the base directory for symlinks
BASE_DIR="data"

# Create an array of symlinks
FILES=($(ls ${BASE_DIR}/*_R1_001.fastq.gz))

# Get the R1 and R2 files for the current task
Sample1=${FILES[$((SLURM_ARRAY_TASK_ID - 1))]}
Sample2=$(echo ${Sample1} | sed 's/_R1_/_R2_/')

# Define the output directory
OUT_DIR="results/2_fastp"
mkdir -p "${OUT_DIR}"

# Extract sample name from the input file path
SAMPLE_NAME=$(basename ${Sample1} | sed 's/_R1_001.fastq.gz//')

# Run fastp with best practices
fastp -i ${Sample1} -I ${Sample2} \
      -o ${OUT_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz \
      -O ${OUT_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz \
      -h ${OUT_DIR}/${SAMPLE_NAME}_fastp_report.html \
      -j ${OUT_DIR}/${SAMPLE_NAME}_fastp_report.json \
      -q 20 -u 30 -n 5 \
      --cut_front --cut_tail --cut_right -w 3

echo -e "\nInfo: fastp job completed for ${SAMPLE_NAME} on $(date) on $(hostname).\n"

