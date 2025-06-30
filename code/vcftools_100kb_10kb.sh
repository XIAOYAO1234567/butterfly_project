#!/bin/bash
#$ -cwd                             
#$ -V                               
#$ -l h_vmem=6G                     
#$ -l h_rt=3:00:00                  
#$ -pe smp 2                       
#$ -N run_piD_allsites_100kb        

# Load vcftools module
module load vcftools

# Define the path to the merged and filtered all-sites VCF file
VCF="bcftools_allsites_filtered_allchr/merged_allchr.Q30_DP16-66_MQ30.vcf.gz"

# Define output directory and prefix
OUTDIR="piD_output_allchr_100kb"
PREFIX="${OUTDIR}/merged_allchr_100kb"

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

echo " Calculating nucleotide diversity (π) using 100 kb non-overlapping windows "
vcftools --gzvcf "$VCF" \
         --window-pi 100000 \
         --window-pi-step 100000 \
         --out "${PREFIX}.pi"

echo " Calculating Tajima's D using 100 kb windows "
vcftools --gzvcf "$VCF" \
         --TajimaD 100000 \
         --out "${PREFIX}.tajima"

echo " Finished: π and Tajima's D results saved in ${OUTDIR}/ "

#!/bin/bash
#$ -cwd                             
#$ -V                               
#$ -l h_vmem=6G                     
#$ -l h_rt=3:00:00                  
#$ -pe smp 2                        
#$ -N run_piD_allsites_10kb         

# Load vcftools module
module load vcftools

# Define the path to the merged and filtered all-sites VCF file
VCF="bcftools_allsites_filtered_allchr/merged_allchr.Q30_DP16-66_MQ30.vcf.gz"

# Define output directory and prefix
OUTDIR="piD_output_allchr_10kb"
PREFIX="${OUTDIR}/merged_allchr_10kb"

# Create output directory if it does not exist
mkdir -p "$OUTDIR"

echo " Calculating nucleotide diversity (π) using 10 kb non-overlapping windows "
vcftools --gzvcf "$VCF" \
         --window-pi 10000 \
         --window-pi-step 10000 \
         --out "${PREFIX}.pi"

echo "=== Calculating Tajima's D using 10 kb windows ==="
vcftools --gzvcf "$VCF" \
         --TajimaD 10000 \
         --out "${PREFIX}.tajima"

echo " Finished: π and Tajima's D results saved in ${OUTDIR}/ "

