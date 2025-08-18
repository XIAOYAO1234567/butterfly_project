#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=4G

start=$(date +%s)

source ~/.bashrc
conda activate rnaseq_qc

mkdir -p multiqc_output
multiqc qc/ -o multiqc_output

end=$(date +%s)
runtime=$((end - start))

echo "MultiQC finished"
echo "Runtime: $runtime seconds"
echo "$(date)"
