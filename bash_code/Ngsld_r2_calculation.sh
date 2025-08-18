#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_vmem=32G
#$ -l h_rt=48:00:00

set -euo pipefail

echo "Job started at $(date)"

module load gcc/12.2.0
module load gsl/2.7.1-gcc-12.2.0
export LD_LIBRARY_PATH=/share/apps/rocky9/spack/apps/linux-rocky9-x86_64_v4/gcc-12.2.0/gsl/2.7.1-ynxvfrc/lib:$LD_LIBRARY_PATH

source ~/.bashrc
conda activate angsd_env

BAM_LIST="bam.filelist"
REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"
REGION="$1"
NGSLD_EXEC="$HOME/software/ngsLD/ngsLD"
NGSLD_OUTDIR="ngsld_output"

SAFE_REGION=$(echo "$REGION" | sed 's/[:\-]/_/g')
OUT_PREFIX="ngsld_${SAFE_REGION}"

mkdir -p "$NGSLD_OUTDIR" logs

echo "Processing region: $REGION"
echo "Output prefix: $OUT_PREFIX"

angsd -bam "$BAM_LIST" -ref "$REF" -r "$REGION" -GL 1 -doGlf 2 -doMajorMinor 1 -doMaf 1 -minMapQ 30 -minQ 30 -P 1 -out "$OUT_PREFIX"

zcat "$OUT_PREFIX.beagle.gz" | tail -n +2 | awk -F"\t" '{split($1,a,"_"); print a[1]"_"a[2]"\t"a[3]}' > "$OUT_PREFIX.pos"
sort -k1,1 -k2,2n "$OUT_PREFIX.pos" -o "$OUT_PREFIX.pos"
sed -i '/^\s*$/d' "$OUT_PREFIX.pos"
sed -i 's/[ \t]*$//' "$OUT_PREFIX.pos"
dos2unix "$OUT_PREFIX.pos" >/dev/null 2>&1 || true
zcat "$OUT_PREFIX.beagle.gz" | tail -n +2 | cut -f2- | gzip > "$OUT_PREFIX.noheader.beagle.gz"

POS_LINES=$(wc -l < "$OUT_PREFIX.pos")
GENO_LINES=$(zcat "$OUT_PREFIX.noheader.beagle.gz" | wc -l)
echo "Number of SNP positions: $POS_LINES"
echo "Number of genotype lines: $GENO_LINES"
if [ "$POS_LINES" -ne "$GENO_LINES" ]; then
  echo "ERROR: SNP positions and genotype lines count mismatch for region $REGION!"
  exit 1
fi

"$NGSLD_EXEC" --geno "$OUT_PREFIX.noheader.beagle.gz" --pos "$OUT_PREFIX.pos" --probs --n_ind 5 --n_sites "$POS_LINES" --max_kb_dist 20 --min_maf 0.05 --out "$NGSLD_OUTDIR/$OUT_PREFIX"

echo "Finished region $REGION at $(date)"
