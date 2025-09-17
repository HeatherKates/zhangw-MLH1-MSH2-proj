#!/bin/bash
#SBATCH --job-name=fastp_fix
#SBATCH --output=logs/fastp_fix_%A_%a.out
#SBATCH --error=logs/fastp_fix_%A_%a.err
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett-b

module load fastp

# Define samples
samples=(
"MLH1KO-2 ../data/MLH1KO-2_S2_L002_R1_001.fastq.gz ../data/MLH1KO-2_S2_L002_R2_001.fastq.gz"
"MLH1KO-3 ../data/MLH1KO-3_S3_L002_R1_001.fastq.gz ../data/MLH1KO-3_S3_L002_R2_001.fastq.gz"
"MLH1R4-3 ../data/MLH1R4-3_S6_L002_R1_001.fastq.gz ../data/MLH1R4-3_S6_L002_R2_001.fastq.gz"
)


# Select sample based on array ID
entry=(${samples[$SLURM_ARRAY_TASK_ID]})
sample_name=${entry[0]}
r1=${entry[1]}
r2=${entry[2]}

# Output files (to fixed_fastqs/)
out1=fixed_fastqs/${sample_name}_R1_fixed.fastq.gz
out2=fixed_fastqs/${sample_name}_R2_fixed.fastq.gz

# Run fastp
fastp \
  --in1 ${r1} \
  --in2 ${r2} \
  --out1 ${out1} \
  --out2 ${out2} \
  --dont_eval_duplication \
  --disable_length_filtering \
  --thread 4

