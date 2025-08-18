#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00
#$ -l h_vmem=32G
#$ -pe smp 1
#$ -N run_chrwise_cov_100kb

source ~/.bashrc
conda activate bedtools_env

CHR_LIST=(
  NC_069083.1 NC_069084.1 NC_069085.1 NC_069086.1 NC_069087.1 NC_069088.1 NC_069089.1 NC_069090.1
  NC_069091.1 NC_069092.1 NC_069093.1 NC_069094.1 NC_069095.1 NC_069096.1 NC_069097.1 NC_069098.1
  NC_069099.1 NC_069100.1 NC_069101.1 NC_069102.1 NC_069103.1 NC_069104.1 NC_069105.1 NC_069106.1
  NC_069107.1 NC_069108.1 NC_069109.1 NC_069110.1
)

BAM_DIR="mapping_results"
OUT_BAM_DIR="mapping_results_chrwise"
WIN_DIR="reference_data/chr_windows"
COV_DIR="coverage_output_chrwise"
REF_IDX="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna.fai"

mkdir -p "$OUT_BAM_DIR" "$WIN_DIR" "$COV_DIR"

for BAM in $BAM_DIR/*.sorted.bam; do
  SAMPLE=$(basename "$BAM" .sorted.bam)
  for CHR in "${CHR_LIST[@]}"; do
    WIN_BED="$WIN_DIR/${CHR}.100kb.bed"
    Q1BAM="$OUT_BAM_DIR/${SAMPLE}.${CHR}.bam"
    COV_FILE="$COV_DIR/${SAMPLE}_${CHR}_coverage.bed"
    if [ ! -f "$WIN_BED" ]; then
      CHR_LEN=$(awk -v chr="$CHR" '$1 == chr {print $2}' "$REF_IDX")
      bedtools makewindows -g <(echo -e "$CHR\t$CHR_LEN") -w 100000 -s 100000 > "$WIN_BED"
    fi
    samtools view -b -q 1 "$BAM" "$CHR" > "$Q1BAM"
    samtools index "$Q1BAM"
    bedtools coverage -a "$WIN_BED" -b "$Q1BAM" -mean > "$COV_FILE"
  done
done

echo "[ALL DONE] Coverage calculation complete (100kb)"
