#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -t 1-28
#$ -N allsites_filtered_allchr

REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"
BAM_DIR="mapping_results"
CHR_LIST="chr.list"
OUT_ROOT="bcftools_allsites_filtered_allchr"

module purge
module load bcftools/1.19-gcc-12.2.0
module load samtools/1.19.2-python-3.11.7-gcc-12.2.0

START=$(date +%s)

# Build chr.list if needed
if [[ ! -f ${CHR_LIST} ]]; then
    echo "[INFO] Creating chr.list from ${REF}.fai"
    samtools faidx "$REF"
    cut -f1 "${REF}.fai" > "${CHR_LIST}"
fi

mapfile -t CHR_ARR < "${CHR_LIST}"
NCHR=${#CHR_ARR[@]}
IDX=$((SGE_TASK_ID - 1))

if (( IDX >= NCHR )); then
    echo "[WARN] Task ${SGE_TASK_ID} exceeds chromosome list length ($NCHR) – exiting"
    exit 0
fi

CHR=${CHR_ARR[$IDX]}
OUT_DIR="${OUT_ROOT}/${CHR}"
VCF_RAW="${OUT_DIR}/merged.${CHR}.allsites.raw.vcf.gz"
VCF_FILT="${OUT_DIR}/merged.${CHR}.allsites.Q30_DP16-66_MQ30.vcf.gz"

mkdir -p "$OUT_DIR"
echo "[INFO] (${SGE_TASK_ID}/${NCHR}) Generating all-sites VCF for ${CHR}"

BAM_LIST=$(ls ${BAM_DIR}/*.sorted.bam 2>/dev/null | tr '\n' ' ')
[[ -z ${BAM_LIST} ]] && { echo "[ERROR] No BAM files found in ${BAM_DIR}"; exit 1; }

# Step 1: mpileup + call (include invariant sites)
bcftools mpileup -Ou -f "$REF" -r "$CHR" ${BAM_LIST} -q 20 -Q 20 -a DP,AD \
  | bcftools call -m -Ou \
  | bcftools view -Oz -o "$VCF_RAW"
bcftools index -f "$VCF_RAW"

# Step 2: apply Q30_DP16-66_MQ30 filtering (keep invariant sites)
bcftools filter \
  -i 'QUAL>=30 && FMT/DP>=16 && FMT/DP<=66 && INFO/MQ>=30' \
  "$VCF_RAW" -Oz -o "$VCF_FILT"
bcftools index -f "$VCF_FILT"

END=$(date +%s)
echo "[DONE] ${CHR}  →  ${VCF_FILT}"
echo "Runtime: $((END - START)) seconds   Completed: $(date)"
