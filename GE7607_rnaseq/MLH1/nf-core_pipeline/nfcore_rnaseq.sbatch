#!/bin/bash
#SBATCH --job-name=MLH1rna.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=120gb
#SBATCH --qos=timgarrett-b
#SBATCH --account=timgarrett
#SBATCH --time=48:00:00
#SBATCH --output=logs/nf-core_rnaseq.resume.out
#SBATCH --error=logs/nf-core_rnaseq.resume.err

module load conda

conda activate rnaseq_env

nextflow run nf-core/rnaseq -c nextflow.config --input fixed_samplesheet.csv \
--outdir OUTPUT --genome GRCh38 --igenomes_ignore \
-profile singularity --mito_name MT \
--fasta /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.112.gtf \
--read_length 150 --save_reference \
--aligner star_rsem --skip_preseq false --trimmer fastp -resume

