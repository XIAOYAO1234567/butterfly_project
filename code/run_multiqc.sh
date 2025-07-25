######## Run MultiQC  ########
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=1:0:0
#$ -l h_vmem=4G

# Starting timestamp
start=$(date +%s)

# Activate Conda environment
source ~/.bashrc
conda activate rnaseq_qc

mkdir -p multiqc_output

# Run MultiQC only on FastQC reports
multiqc qc/ -o multiqc_output

end=$(date +%s)
runtime=$((end - start))

echo "MultiQC (FastQC-only) finished"
echo "Runtime: $runtime seconds"
echo $(date)

