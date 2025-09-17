#!/bin/bash
set -euo pipefail

# Define target directory
TARGET_DIR="/orange/zhangw/GENOME_REFS/Mus_musculus/GenomicVariants/GRCm38"

# Create directory if it doesn't exist
mkdir -p "$TARGET_DIR"
cd "$TARGET_DIR"

# List of URLs to download
URLS=(
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz.tbi"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz"
  "https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.indels.dbSNP142.normed.vcf.gz.tbi"
)

# Download each file
for url in "${URLS[@]}"; do
  echo "Downloading $url..."
  wget -N "$url"
done

echo "All files downloaded to $TARGET_DIR"
