#!/bin/bash
#SBATCH --job-name=makeindex
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=16gb
#SBATCH --qos=zhangw-b
#SBATCH --account=zhangw
#SBATCH -t 04:00:00
#SBATCH --output=logs/makeindex.out
#SBATCH --error=logs/makeindex.err

# Load Salmon module
module load salmon/1.10.1

# Define the directory for the Salmon index
path="../salmon_index"

# Check if the Salmon index directory exists
if [ -d "$path" ]; then
    echo "The directory $path already exists. Exiting without executing further commands."
    exit 0
fi

# Create the directory for the Salmon index
mkdir -p $path

# Download the transcriptome and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip $path/Homo_sapiens.GRCh38.cdna.all.fa.gz

# Download the GTF file and move it to the designated directory
wget -P $path ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
gunzip $path/Homo_sapiens.GRCh38.104.gtf.gz

# Create the Salmon index in the designated directory
salmon index -t $path/Homo_sapiens.GRCh38.cdna.all.fa -i $path/salmon_index_human
