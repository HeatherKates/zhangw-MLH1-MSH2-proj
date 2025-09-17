#!/bin/bash
#SBATCH --job-name=MSH2rna.%A_%a
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=120gb
#SBATCH --qos=cancercenter-dept-b
#SBATCH --account=cancercenter-dept
#SBATCH --time=48:00:00
#SBATCH --output=logs/nf-core_rnaseq.skiptrim.out
#SBATCH --error=logs/nf-core_rnaseq.skiptrim.err

module load conda

conda activate rnaseq_env

nextflow run nf-core/rnaseq -c nextflow.config --input fixed_samplesheet.csv \
--outdir OUTPUT --genome GRCh38 --igenomes_ignore \
-profile singularity --mito_name MT \
--fasta /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--gtf /orange/cancercenter-dept/GENOMES/iGenomes/references/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/Homo_sapiens.GRCh38.112.gtf \
--read_length 150 --save_reference \
--aligner star_rsem --skip_preseq false --skip_trimming true -resume --save_bam_intermediates true --bam_csi_index true


