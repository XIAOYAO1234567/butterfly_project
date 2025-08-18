#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 4
#$ -l h_rt=12:00:00
#$ -l h_vmem=8G
#$ -t 1-28

REF="reference_data/GCF_947172395.1_ilBicAnyn1.1_genomic.fna"
BAM_DIR="mapping_results"
CHR_LIST="chr.list"
OUT_ROOT="bcftools_ml_allchr_10-70"

module purge
module load bcftools/1.19-gcc-12.2.0
module load samtools/1.19.2-python-3.11.7-gcc-12.2.0

START=$(date +%s)

if [[ ! -f ${CHR_LIST} ]]; then
    samtools faidx "$REF"
    cut -f1 "${REF}.fai" > "${CHR_LIST}"
fi

mapfile -t CHR_ARR < "${CHR_LIST}"
NCHR=${#CHR_ARR[@]}
IDX=$((SGE_TASK_ID - 1))

if (( IDX >= NCHR )); then
    exit 0
fi

CHR=${CHR_ARR[$IDX]}
OUT_DIR="${OUT_ROOT}/${CHR}"
VCF_RAW="${OUT_DIR}/merged.${CHR}.snps.raw.vcf.gz"
VCF_FILT="${OUT_DIR}/merged.${CHR}.snps.Q30_DP10-70_MQ30.vcf.gz"

mkdir -p "$OUT_DIR"

BAM_LIST=$(ls ${BAM_DIR}/*.sorted.bam 2>/dev/null | tr '\n' ' ')
[[ -z ${BAM_LIST} ]] && exit 1

bcftools mpileup -Ou -f "$REF" -r "$CHR" ${BAM_LIST} -q 20 -Q 20 -a DP,AD \
  | bcftools call -mv -Oz -o "$VCF_RAW"
bcftools index -f "$VCF_RAW"

bcftools filter \
  -i 'QUAL>=30 && FMT/DP>=10 && FMT/DP<=70 && INFO/MQ>=30' \
  "$VCF_RAW" -Oz -o "$VCF_FILT"
bcftools index -f "$VCF_FILT"

END=$(date +%s)
echo "[DONE] ${CHR} → ${VCF_FILT}"
echo "Runtime: $((END - START)) seconds   Completed: $(date)"
