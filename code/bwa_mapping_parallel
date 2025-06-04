######## Mapping with BWA ########
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=240:0:0
#$ -l h_vmem=12G
#$ -N run_bwa_single

# Starting timestamp
start=$(date +%s)

# Activate Conda environment with bwa & samtools
source ~/.bashrc
conda activate genome_mapping

# Input sample name passed from qsub
SAMPLE=$1

# Set reference genome
REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"

# Input files
FQ1="clean_data/${SAMPLE}_1.clean.fastq.gz"
FQ2="clean_data/${SAMPLE}_2.clean.fastq.gz"

# Output files
SAM="mapping_results/${SAMPLE}.sam"
BAM="mapping_results/${SAMPLE}.bam"
SORTED_BAM="mapping_results/${SAMPLE}.sorted.bam"

# Index the reference if not done already
if [ ! -f "${REF}.bwt" ]; then
  echo "[INFO] Indexing reference genome..."
  bwa index "$REF"
fi

# Run BWA MEM
bwa mem -t 4 "$REF" "$FQ1" "$FQ2" > "$SAM"

# Convert to BAM
samtools view -@ 4 -Sb "$SAM" > "$BAM"
rm "$SAM"

# Sort BAM
samtools sort -@ 4 "$BAM" -o "$SORTED_BAM"
rm "$BAM"

# Index BAM
samtools index "$SORTED_BAM"

# Finishing timestamp
end=$(date +%s)
runtime=$((end - start))

echo "Finished mapping $SAMPLE"
echo "Runtime: $runtime seconds"
echo $(date)




#!/bin/bash
for fq1 in clean_data/*_1.clean.fastq.gz; do
  sample=$(basename "$fq1" _1.clean.fastq.gz)
  qsub run_bwa_single_sample.sh "$sample"
done
