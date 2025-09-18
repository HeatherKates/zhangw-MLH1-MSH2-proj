#!/bin/bash
#SBATCH --job-name=rnavar
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3gb
#SBATCH --qos=cancercenter-dept-b
#SBATCH --account=cancercenter-dept
#SBATCH -t 48:00:00
#SBATCH --output=logs/rnavar.%j.out
#SBATCH --error=logs/rnavar.%j.err

module load nextflow
module load singularity

# Optional: store singularity containers in a clean cache
export NXF_SINGULARITY_CACHEDIR=$PWD/singularity_cache

nextflow run nf-core/rnavar --outdir OUTPUT  --input samplesheet.csv \
  --genome GRCm38 \
  -profile singularity \
  -c nextflow.config \
  --dbsnp /orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38/mgp.v5.merged.snps_all.dbSNP142.vcf.gz \
  --dbsnp_tbi /orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi \
  --known_indels /orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz \ 
  --known_indels_tbi /orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz.tbi\
  -resume
