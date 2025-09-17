#!/bin/bash
#SBATCH --job-name=salmon_pipeline
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16gb
#SBATCH --qos=REPLACE_WITH_YOUR_QOS
#SBATCH --account=REPLACE_WITH_YOUR_ACCOUNT
#SBATCH -t 08:00:00
#SBATCH --output=logs/salmon.%j.out
#SBATCH --error=logs/salmon.%j.err
#SBATCH --array=1-12

echo "This is task $SLURM_ARRAY_TASK_ID"

echo -e "\nInfo: Starting a job on $(date) on $(hostname) in $(pwd).\n"

###################################################################
# This file is to perform transcript quantification with Salmon
###################################################################

module load salmon/1.10.1

# Define the base directory
BASE_DIR="results/2_fastp"

# Create an array of trimmed files dynamically
FILES=$(ls ${BASE_DIR}/*_R1_trimmed.fastq.gz)

# Get the files for the current task
FILE=$(echo $FILES | tr ' ' '\n' | sed -n "${SLURM_ARRAY_TASK_ID}p")

# Define the R1 and R2 files
Sample1=${FILE}
Sample2=$(echo ${FILE} | sed 's/_R1_/_R2_/')

# Define the output directory
OUT_DIR="results/salmon"

# Create output directory if it doesn't exist
mkdir -p ${OUT_DIR}

# Extract sample name from the input file path
SAMPLE_NAME=$(basename ${Sample1} | sed 's/_R1_trimmed.fastq.gz//')

# Path to the Salmon index
SALMON_INDEX="../salmon_index/salmon_index_human"

# Run Salmon quantification with quasi-mapping
salmon quant -i ${SALMON_INDEX} -l A -1 ${Sample1} -2 ${Sample2} -p 4 -o ${OUT_DIR}/${SAMPLE_NAME}_quant --seqBias --gcBias

echo -e "\nInfo: Salmon quantification job completed on $(date) on $(hostname).\n"

