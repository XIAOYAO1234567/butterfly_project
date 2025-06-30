#!/bin/bash
#$ -N run_angsd_theta
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_vmem=32G
#$ -l h_rt=12:00:00

set -euo pipefail
echo "Starting ANGSD theta analysis at $(date)"

# Activate conda environment
source ~/.bashrc
conda activate angsd_env

# Optional: show modules
module list || true

# Paths & Parameters 
BAM_LIST="bam.filelist"
REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"
CHR_INCLUDE="chr.list"
OUTDIR="angsd_output"
PREFIX="theta"

mkdir -p "$OUTDIR"

echo "[INFO] Removing previous outputs if exist"
rm -f "$OUTDIR/$PREFIX."*

# Step 1: Compute SAF 
echo "[Step 1] Running ANGSD to compute SAF"
angsd -bam "$BAM_LIST" \
  -ref "$REF" -anc "$REF" \
  -rf "$CHR_INCLUDE" \
  -out "$OUTDIR/$PREFIX" \
  -doSaf 1 \
  -GL 1 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doCounts 1 \
  -minMapQ 30 -minQ 30 \
  -setMinDepthInd 16 -setMaxDepthInd 66 \
  -minInd 3 \
  -P 4

# Step 2: Estimate folded SFS
echo "[Step 2] Estimating folded SFS with realSFS"
realSFS "$OUTDIR/${PREFIX}.saf.idx" -P 4 -fold 1 > "$OUTDIR/${PREFIX}.sfs"

# Step 3: Compute thetas 
echo "[Step 3] Computing theta with saf2theta"
realSFS saf2theta "$OUTDIR/${PREFIX}.saf.idx" \
  -sfs "$OUTDIR/${PREFIX}.sfs" \
  -outname "$OUTDIR/$PREFIX"

# Step 4: Print site-wise theta values 
echo "[Step 4a] Printing site-wise theta values"
thetaStat print "$OUTDIR/${PREFIX}.thetas.idx" > "$OUTDIR/${PREFIX}.thetas.print.txt"

# Step 5: Sliding window analysis
echo "[Step 4b] Calculating sliding window stats"
thetaStat do_stat "$OUTDIR/${PREFIX}.thetas.idx" -win 10000 -step 10000 -outnames "$OUTDIR/${PREFIX}.10kb"
thetaStat do_stat "$OUTDIR/${PREFIX}.thetas.idx" -win 100000 -step 100000 -outnames "$OUTDIR/${PREFIX}.100kb"

echo "ANGSD theta analysis completed successfully at $(date)"
