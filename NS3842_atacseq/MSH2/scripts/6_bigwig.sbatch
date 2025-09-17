#!/bin/bash
#SBATCH --job-name=bigwig.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=100gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=02:00:00
#SBATCH --output=../logs/bigwig_%A_%a.out
#SBATCH --error=../logs/bigwig_%A_%a.err
#SBATCH --array=0-5

module load deeptools/3.5.2
module load samtools/1.20

# Get list of BAM files
bam_files=(results/4a_samtools/*.sorted.bam)

# Select the BAM file for the current array task
bam_file=${bam_files[$SLURM_ARRAY_TASK_ID]}

# Extract the sample name from the BAM file path (remove directory and file extension)
sample_name=$(basename $bam_file .sorted.bam)

# Sort the filtered BAM by coordinate (if necessary) for visualization and diffbind
samtools sort -o results/6_bigwig/${sample_name}_sort_n.bam $bam_file

# Re-index the BAM file just in case it's needed
samtools index results/6_bigwig/${sample_name}_sort_n.bam

# Create the BigWig file using bamCoverage
bamCoverage --bam results/6_bigwig/${sample_name}_sort_n.bam -o results/6_bigwig/${sample_name}.bw --binSize 10 --normalizeUsing RPKM -p=6

