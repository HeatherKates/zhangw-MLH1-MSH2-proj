#!/bin/bash

# Define the input file containing the file paths
input_file="file_paths.txt"
# Define the output CSV file
output_file="file_metadata.csv"

# Write the header to the output file
echo "File,Analysis,SampleType,FileType,Software,Description" > $output_file

# Read each file path and process it
while read -r filepath; do
    # Extract file name
    filename=$(basename "$filepath")

    # Determine analysis type based on the path
    if [[ "$filepath" == *"ATACseq"* ]]; then
        analysis="ATACseq"
    elif [[ "$filepath" == *"CUTANDRUN"* ]]; then
        analysis="CUT&RUN"
    elif [[ "$filepath" == *"RNAseq"* ]]; then
        analysis="RNAseq"
    else
        analysis="Unknown"
    fi

    # Determine sample type based on directory structure
    if [[ "$filepath" == *"MLH1"* ]]; then
        sample_type="MLH1"
    elif [[ "$filepath" == *"MSH2"* ]]; then
        sample_type="MSH2"
    else
        sample_type="Unknown"
    fi

    # Guess file type based on extension
    case "$filename" in
        *.bam) file_type="BAM" ;;
        *.bai) file_type="BAI" ;;
        *.csv) file_type="CSV" ;;
        *.xlsx) file_type="Excel" ;;
        *.bed) file_type="BED" ;;
        *.bedgraph) file_type="BedGraph" ;;
        *.RDATA) file_type="RData" ;;
        *.narrowPeak) file_type="NarrowPeak" ;;
        *) file_type="Other" ;;
    esac

    # Attempt to infer software based on directory or file type
    if [[ "$filepath" == *"bigwig"* ]]; then
        software="BigWig Tools"
    elif [[ "$filepath" == *"samtools"* ]]; then
        software="Samtools"
    elif [[ "$filepath" == *"macs2"* ]]; then
        software="MACS2"
    elif [[ "$filepath" == *"seacr"* ]]; then
        software="SEACR"
    elif [[ "$filepath" == *"deseq2"* ]]; then
        software="DESeq2"
    else
        software="Unknown"
    fi

    # Generate a basic description
    description="${file_type} file generated during ${analysis} analysis."

    # Append the line to the output CSV file
    echo "$filepath,$analysis,$sample_type,$file_type,$software,$description" >> $output_file

done < "$input_file"

# Notify completion
echo "Metadata extraction complete. Results saved to $output_file."

