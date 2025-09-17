#!/bin/bash
#SBATCH --job-name=MSH2_atacseq
#SBATCH --output=logs/atacseq_%j.out
#SBATCH --error=logs/atacseq_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=240G
#SBATCH --account=timgarrett
#SBATCH --qos=timgarrett-b

# Load required modules
module load nextflow/24.10.5
module load apptainer/1.3.5
module load java/11.0.1

# Run nf-core/atacseq pipeline

nextflow run nf-core/atacseq -c nextflow.config --input samplesheet.csv \
--outdir OUTPUT --genome GRCh38 --igenomes_ignore \
--blacklist /orange/cancercenter-dept/resources/BLACKLISTS/hg38-blacklist.v2.ensembl.bed \
-profile singularity --mito_name MT \
--fasta /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.112.gtf \
--save_macs_pileup --aligner bowtie2 --read_length 150 --save_reference

