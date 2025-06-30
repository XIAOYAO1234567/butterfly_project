######## Adapter & Quality Trimming ########
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=240:0:0
#$ -l h_vmem=8G

# Starting timestamp
start=$(date +%s)

# Activate Conda environment
source ~/.bashrc
conda activate rnaseq_qc

# Create output directories
mkdir -p clean_data fastp_reports

# Run fastp for each paired-end sample
for sample in raw_data/*_1.fastq.gz
do
  base=$(basename "$sample" _1.fastq.gz)

  fastp \
    -i raw_data/${base}_1.fastq.gz \
    -I raw_data/${base}_2.fastq.gz \
    -o clean_data/${base}_1.clean.fastq.gz \
    -O clean_data/${base}_2.clean.fastq.gz \
    --thread 4 \
    --detect_adapter_for_pe \
    --trim_front1 10 --trim_tail1 10 \
    --trim_front2 10 --trim_tail2 10 \
    --length_required 50 \
    --html fastp_reports/${base}.html \
    --json fastp_reports/${base}.json \
    --report_title "${base} fastp report"
done

# Finishing timestamp
end=$(date +%s)
runtime=$((end - start))

echo Finished
echo "Runtime: $runtime seconds"
echo $(date)
