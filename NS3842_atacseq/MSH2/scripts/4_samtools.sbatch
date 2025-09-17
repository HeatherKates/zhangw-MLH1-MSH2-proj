#!/bin/bash
#SBATCH --job-name=samtools.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=24:00:00
#SBATCH --output=../logs/samtools_%A_%a.out
#SBATCH --error=../logs/samtools_%A_%a.err
#SBATCH --array=0-5  # Assuming 6 samples in the list

# Load the samtools module
module load samtools/1.20

# Read the sample name from the sample list
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" MSH2_sample_list.txt)

# Step 1: Convert the SAM file to BAM
samtools view -h -b results/3_bowtie/${sample}.sam > results/4a_samtools/${sample}.bam

# Step 2: Apply filters to keep only mapped, properly paired reads with high quality, and save to temp BAM
samtools view -h -b -F 1804 -f 2 -q 30 -@ 4 results/4a_samtools/${sample}.bam > results/4a_samtools/${sample}.filtered.bam

# Step 3: Use awk to exclude non-standard chromosomes and output a filtered BAM
samtools view -h results/4a_samtools/${sample}.filtered.bam | \
    awk '$1 ~ /^@/ || $3 ~ /^[1-9XY][0-9]*$/' | \
    samtools view -b - > results/4a_samtools/${sample}.filt.bam

# Step 4: Sort the filtered BAM by read name (used for peak calling)
samtools sort -n -O BAM -@ 4 results/4a_samtools/${sample}.filt.bam > results/4a_samtools/${sample}.sorted.bam

# Step 5: Run flagstat on aligned BAM to get alignment statistics
mkdir -p results/4b_flagstat
samtools flagstat -@ 4 results/4a_samtools/${sample}.bam > results/4b_flagstat/${sample}.bam.flagstat.log

