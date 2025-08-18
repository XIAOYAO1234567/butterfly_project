#!/bin/bash
set -euo pipefail

VCF="merged.NC_069094.1.snps.Q30_DP10-70_MQ30.vcf.gz" # we can change the chromosome in this place, 094 chromosome is just an example
CHR="NC_069094.1" 

BAL_REGIONS=("370001 400000" "680001 710000") # we can change the regions in this place, it is just an example here
NEU_REGIONS=("530001 580000")

rm -f balancing.body.vcf neutral.body.vcf header.vcf

for reg in "${BAL_REGIONS[@]}"; do
  read -r start end <<< "$reg"
  zcat "$VCF" | awk -v chr="$CHR" -v s="$start" -v e="$end" '$1==chr && $2>=s && $2<=e' >> balancing.body.vcf
done

for reg in "${NEU_REGIONS[@]}"; do
  read -r start end <<< "$reg"
  zcat "$VCF" | awk -v chr="$CHR" -v s="$start" -v e="$end" '$1==chr && $2>=s && $2<=e' >> neutral.body.vcf
done

sed -i 's/DP=/MT=3;DP=/' balancing.body.vcf
sed -i 's/DP=/MT=1;DP=/' neutral.body.vcf

zcat "$VCF" | awk '/^##/ || /^#CHROM/' > header.vcf

cat header.vcf neutral.body.vcf > NC_069094.1_MT1_neutral.vcf
cat header.vcf balancing.body.vcf > NC_069094.1_MT3_balancing.vcf

for input in NC_069094.1_MT*.vcf; do
  output="${input%.vcf}.cleaned.vcf"
  awk 'BEGIN{OFS="\t"}
       { if($0~/^#/) {print; next}
         if($4!~/^[ACGT]$/ || $5!~/^[ACGT]$/) next
         if($8~/AC=\./ || $8~/AN=\./) next
         for(i=10;i<=NF;i++){ gsub(/\.\/\.|\.\|\.|\./,"0/0",$i) }
         print }' "$input" > "$output"
done
