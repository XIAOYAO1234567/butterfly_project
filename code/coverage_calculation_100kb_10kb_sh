#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=24:00:00
#$ -l h_vmem=32G
#$ -pe smp 1
#$ -N run_chrwise_cov_all

# Load conda environment with samtools and bedtools
source ~/.bashrc
conda activate bedtools_env

# Define variables
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

# Loop through all samples 
for BAM in $BAM_DIR/*.sorted.bam; do
  SAMPLE=$(basename "$BAM" .sorted.bam)
  echo "[INFO] Processing sample: $SAMPLE"

  for CHR in "${CHR_LIST[@]}"; do
    echo "  Processing: $CHR"

    WIN_BED="$WIN_DIR/${CHR}.100kb.bed"
    Q1BAM="$OUT_BAM_DIR/${SAMPLE}.${CHR}.bam"
    COV_FILE="$COV_DIR/${SAMPLE}_${CHR}_coverage.bed"

    # Generate window file if not exist
    if [ ! -f "$WIN_BED" ]; then
      CHR_LEN=$(awk -v chr="$CHR" '$1 == chr {print $2}' "$REF_IDX")
      if [ -z "$CHR_LEN" ]; then
        echo "[WARNING] Chromosome $CHR not found in index. Skipping."
        continue
      fi
      bedtools makewindows -g <(echo -e "$CHR\t$CHR_LEN") -w 100000 -s 100000 > "$WIN_BED"
    fi

    # Extract and filter
    samtools view -b -q 1 "$BAM" "$CHR" > "$Q1BAM"
    samtools index "$Q1BAM"

    # Calculate coverage
    bedtools coverage -a "$WIN_BED" -b "$Q1BAM" -mean > "$COV_FILE"

    if [ ! -s "$COV_FILE" ]; then
      echo "[WARNING] Coverage file is empty: $COV_FILE"
    else
      echo "[DONE] $SAMPLE $CHR coverage done."
    fi
  done
done

echo "[ALL DONE] Coverage calculation complete for all chromosomes and samples."




#!/bin/bash
#$ -cwd
#$ -j y
#$ -l h_rt=12:00:00
#$ -l h_vmem=16G
#$ -pe smp 1

# Activate conda environment 
source ~/.bashrc
conda activate bedtools_env

# Define variables
CHR_LIST=(
  NC_069083.1 NC_069084.1 NC_069085.1 NC_069086.1 NC_069087.1 NC_069088.1 NC_069089.1 NC_069090.1
  NC_069091.1 NC_069092.1 NC_069093.1 NC_069094.1 NC_069095.1 NC_069096.1 NC_069097.1 NC_069098.1
  NC_069099.1 NC_069100.1 NC_069101.1 NC_069102.1 NC_069103.1 NC_069104.1 NC_069105.1 NC_069106.1
  NC_069107.1 NC_069108.1 NC_069109.1 NC_069110.1
)

BAM_DIR="mapping_results_chrwise/mapping_results_chrwise"             
WIN_DIR="reference_data/chr_windows_10kb"
COV_DIR="coverage_output_chrwise_10kb"
mkdir -p "$COV_DIR"

# Loop through all samples 
for BAM in ${BAM_DIR}/*.bam; do
  SAMPLE_CHR=$(basename "$BAM" .bam)   # e.g., SRR17425635.NC_069083.1
  SAMPLE=$(echo $SAMPLE_CHR | cut -d. -f1)
  CHR=$(echo $SAMPLE_CHR | cut -d. -f2-3)

  WIN_BED="$WIN_DIR/${CHR}.10kb.bed"
  COV_FILE="$COV_DIR/${SAMPLE}_${CHR}_coverage.bed"

  if [ ! -f "$WIN_BED" ]; then
    echo "[WARNING] Missing window file: $WIN_BED"
    continue
  fi

  echo "[INFO] Calculating 10kb coverage: $SAMPLE - $CHR"
  bedtools coverage -a "$WIN_BED" -b "$BAM" -mean > "$COV_FILE"

  if [ ! -s "$COV_FILE" ]; then
    echo "[WARNING] Empty output: $COV_FILE"
  fi
done
echo "[ALL DONE] 10kb coverage calculation finished!"
