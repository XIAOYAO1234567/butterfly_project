#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 2
#$ -l h_rt=6:00:00
#$ -l h_vmem=4G
#$ -N run_fastp_reports

# Activate conda
source ~/.bashrc
conda activate rnaseq_qc

mkdir -p clean_data_fastp_reports

# generate reports
for fq1 in clean_data/*_1.clean.fastq.gz
do
  sample=$(basename "$fq1" _1.clean.fastq.gz)
  fq2="clean_data/${sample}_2.clean.fastq.gz"

  fastp \
    -i "$fq1" \
    -I "$fq2" \
    --thread 2 \
    --disable_adapter_trimming \
    --disable_quality_filtering \
    --disable_length_filtering \
    --html "clean_data_fastp_reports/${sample}_clean.html" \
    --json "clean_data_fastp_reports/${sample}_clean.json" \
    --report_title "${sample} Clean Fastp Report"
done
