#!/bin/bash
#SBATCH --job-name=trim_galore.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=20gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --output=../logs/trim_galore_%A_%a.out
#SBATCH --error=../logs/trim_galore_%A_%a.err
#SBATCH --array=0-5 # For 6 samples in the list

# Load required modules
module load trim_galore
module load fastqc

# Define the directory containing the symlinks to individual fastq.gz files
symlink_dir="data"

# Read the sample name from the sample list
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" MLH1_sample_list.txt)

# Define file paths for R1 and R2 for both lanes
R1_L001="${symlink_dir}/${sample}_L001_R1_001.fastq.gz"
R2_L001="${symlink_dir}/${sample}_L001_R2_001.fastq.gz"
R1_L002="${symlink_dir}/${sample}_L002_R1_001.fastq.gz"
R2_L002="${symlink_dir}/${sample}_L002_R2_001.fastq.gz"

# Create output directories if they don't exist
mkdir -p results/1_trimgalore/
mkdir -p results/2_fastQC/

echo "Running trimgalore on ${R1_L001} ${R2_L001}"
# Run trim_galore on both lanes for R1 and R2
trim_galore --phred33 --paired --fastqc --cores 2 --output_dir results/1_trimgalore/ --basename ${sample}_L001 $R1_L001 $R2_L001

echo "Running trimgalore on  ${R1_L002} ${R2_L002}"
trim_galore --phred33 --paired --fastqc --cores 2 --output_dir results/1_trimgalore/ --basename ${sample}_L002 $R1_L002 $R2_L002

# Run fastqc on the raw input files (R1 and R2 for both lanes)
fastqc -o results/2_fastQC/ $R1_L001 $R2_L001
fastqc -o results/2_fastQC/ $R1_L002 $R2_L002

# Move fastqc outputs into a specific directory
mv results/1_trimgalore/*_fastqc* results/2_fastQC/

