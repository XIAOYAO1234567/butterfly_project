step 1：Define regions and extract VCF subsets

# sweep regions（MT=2）
zcat merged.NC_069085.1.snps.Q30_DP10-70_MQ30.vcf.gz | \
awk -v chr="NC_069085.1" -v start=10725001 -v end=10750000 \
'$1 == chr && $2 >= start && $2 <= end' > sweep_10725001.txt

zcat merged.NC_069085.1.snps.Q30_DP10-70_MQ30.vcf.gz | \
awk -v chr="NC_069094.1" -v start=11040001 -v end=11060000 \
'$1 == chr && $2 >= start && $2 <= end' > sweep_11040001.txt

# neutral regions（MT=1）
zcat merged.NC_069094.1.snps.Q30_DP10-70_MQ30.vcf.gz | \
awk -v chr="NC_069094.1" -v start=10860001 -v end=10920000 \
'$1 == chr && $2 >= start && $2 <= end' > neutral_10860001.txt

step 2：Add MT information for each type of file.

sed -i 's/DP=/MT=2;DP=/g' sweep_10725001.txt
sed -i 's/DP=/MT=2;DP=/g' sweep_11040001.txt
sed -i 's/DP=/MT=1;DP=/g' neutral_10860001.txt

step 3：Concatenate the header and the data body to form a complete VCF.

zcat merged.NC_069094.1.snps.Q30_DP10-70_MQ30.vcf.gz | awk '/^##fileformat/ || /^##INFO/ || /^#CHROM/' > header.vcf

cat header.vcf neutral_10860001.txt > NC_069085.1_MT1_neutral.vcf
cat header.vcf sweep_10725001.txt sweep_11040001.txt > NC_069085.1_MT2_sweep.vcf

step 4:

#!/bin/bash

for input in NC_069085.1_MT*.vcf; do
    output="${input%.vcf}.cleaned.vcf"

    awk '
    BEGIN { OFS = "\t" }
    {
        if ($0 ~ /^#/) { print; next }

        if ($4 !~ /^[ACGT]$/ || $5 !~ /^[ACGT]$/) next
        if ($8 ~ /AC=\./ || $8 ~ /AN=\./) next

        for (i = 10; i <= NF; i++) {
            gsub(/\.\/\.|\.\|\.|\./, "0/0", $i)
        }

        print
    }' "$input" > "$output"
done
