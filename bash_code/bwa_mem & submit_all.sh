#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=240:0:0
#$ -l h_vmem=12G

start=$(date +%s)

source ~/.bashrc
conda activate genome_mapping

SAMPLE=$1
REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"
FQ1="clean_data/${SAMPLE}_1.clean.fastq.gz"
FQ2="clean_data/${SAMPLE}_2.clean.fastq.gz"
SAM="mapping_results/${SAMPLE}.sam"
BAM="mapping_results/${SAMPLE}.bam"
SORTED_BAM="mapping_results/${SAMPLE}.sorted.bam"

if [ ! -f "${REF}.bwt" ]; then
  bwa index "$REF"
fi

bwa mem -t 4 "$REF" "$FQ1" "$FQ2" > "$SAM"
samtools view -@ 4 -Sb "$SAM" > "$BAM"
rm "$SAM"
samtools sort -@ 4 "$BAM" -o "$SORTED_BAM"
rm "$BAM"
samtools index "$SORTED_BAM"

end=$(date +%s)
runtime=$((end - start))

echo "Finished mapping $SAMPLE"
echo "Runtime: $runtime seconds"
echo $(date)


#!/bin/bash
for fq1 in clean_data/*_1.clean.fastq.gz; do
  sample=$(basename "$fq1" _1.clean.fastq.gz)
  qsub run_bwa_single.sh "$sample"
done
