#!/bin/bash

# Set file paths and sample list
base_dir="."  
sample_list="$base_dir/MSH2_sample_list.txt"
output_summary="$base_dir/results/9_summary/ATACseq_results_summary.csv"

# Write header to summary CSV
echo "Sample,Lane,Read,Raw_Reads,Trimmed_Reads,Adapter_Content,GC_Content,Duplicate_Level,Total_Reads_FastQC,Mapped_Reads,Properly_Paired_Reads,Peaks_Called,Differential_Peaks,BigWig_File" > $output_summary

# Loop through each sample
while read sample; do

    # Loop through each lane (L001 and L002) and read (R1 and R2)
    for lane in 001 002; do
        for read in R1 R2; do

            # Initialize variables for each sample, lane, and read
            raw_reads=""
            trimmed_reads=""
            adapter_content=""
            gc_content=""
            duplicate_level=""
            total_reads_fastqc=""
            mapped_reads=""
            properly_paired_reads=""
            peaks_called=""
            differential_peaks=""
            bigwig_file=""

trim_report=$(ls "$base_dir/results/1_trimgalore"/*${sample}*${lane}_${read}_001.fastq.gz_trimming_report.txt 2>/dev/null)
if [ -f "$trim_report" ]; then
	raw_reads=$(grep "Total reads processed:" "$trim_report" | awk '{print $4}'|sed 's/,//g')
	trimmed_reads=$(grep "Reads written (passing filters)" "$trim_report" | awk '{print $5}' | sed 's/,//g')
else
    echo "WARNING: Trimming report not found for $sample $lane $read, expected file: $trim_report" >> "$base_dir/QC_warnings.log"
fi

# Step 2: Find and unzip the FastQC summary file
read_num=$(echo ${read} | cut -c 2)
fastqc_zip=$(ls "$base_dir/results/2_fastQC"/*${sample}*${lane}_val_${read_num}_fastqc.zip 2>/dev/null)
if [ -f "$fastqc_zip" ]; then
    unzip -p "$fastqc_zip" "*fastqc_data.txt" > "${sample}_fastqc_data.txt"

# Parse FastQC data for the correct values
total_reads_fastqc=$(grep "Total Sequences" "${sample}_fastqc_data.txt" | awk '{print $3}'|sed 's/%//g')
adapter_content=$(grep "Adapter Content" "${sample}_fastqc_data.txt" | awk '{print $3}' | sed 's/%//g')
gc_content=$(grep "Per sequence GC content" "${sample}_fastqc_data.txt"| awk '{print $5}')
duplicate_level=$(grep "Total Deduplicated Percentage" "${sample}_fastqc_data.txt" | awk '{print $4}' | sed 's/%//g')
    
    # Remove the extracted file after processing
    rm "${sample}_fastqc_data.txt"
else
    echo "WARNING: FastQC summary not found for $sample $lane $read, looking for file: $fastqc_zip" >> "$base_dir/QC_warnings.log"
fi

            # Step 3: Parse flagstat results for read mapping metrics (constant per sample)
            flagstat_file="$base_dir/results/4b_flagstat/${sample}.bam.flagstat.log"
            if [ -f "$flagstat_file" ]; then
                mapped_reads=$(grep "mapped (" "$flagstat_file" | awk '{print $1}'|head -1)
                properly_paired_reads=$(grep "properly paired (" "$flagstat_file" | awk '{print $1}')
            else
                echo "WARNING: Flagstat file not found for $sample, looking for file: $flagstat_file" >> "$base_dir/QC_warnings.log"
            fi

            # Step 4: Count the number of peaks called by Genrich (constant per sample)
            genrich_peaks="$base_dir/results/5_genrich/${sample}.narrowPeak"
            if [ -f "$genrich_peaks" ]; then
                peaks_called=$(wc -l "$genrich_peaks" | awk '{print $1}')
            else
                echo "WARNING: Genrich peaks file not found for $sample, looking for file: $genrich_peaks" >> "$base_dir/QC_warnings.log"
            fi

            # Step 5: Check for BigWig file existence (constant per sample)
            bigwig_file="$base_dir/results/6_bigwig/${sample}.bw"
            if [ -f "$bigwig_file" ]; then
                bigwig_file=${bigwig_file}
            else
                bigwig_file="No file"
                echo "WARNING: BigWig file not found for $sample, looking for file: $bigwig_file" >> "$base_dir/QC_warnings.log"
            fi

            # Step 6: Check the number of differential peaks from DiffBind results (constant per sample)
# Determine which differential binding results file to use based on sample name
if [[ $sample == *"KO"* ]]; then
    diffbind_results=$(ls "$base_dir"/results/7_diffbind/*KO_filtered_differential_binding_results.csv)
elif [[ $sample == *"R4"* ]]; then
    diffbind_results=$(ls "$base_dir"/results/7_diffbind/*R4_filtered_differential_binding_results.csv)
else
    diffbind_results=""
fi

# Check if the correct file exists and count the differential peaks
if [ -f "$diffbind_results" ]; then
    differential_peaks=$(awk 'END{print NR-1}' "$diffbind_results")  # Exclude header row
else
    echo "WARNING: DiffBind results not found for $sample, looking for file: $diffbind_results" >> "$base_dir/QC_warnings.log"
    differential_peaks="NA"
fi


            # Step 7: Output results for the sample, lane, and read to the summary CSV
            echo "$sample,$lane,$read,$raw_reads,$trimmed_reads,$adapter_content,$gc_content,$duplicate_level,$total_reads_fastqc,$mapped_reads,$properly_paired_reads,$peaks_called,$differential_peaks,$bigwig_file" >> "$output_summary"

        done  # End of read loop (R1, R2)
    done  # End of lane loop (L001, L002)

done < "$sample_list"  # End of sample loop

