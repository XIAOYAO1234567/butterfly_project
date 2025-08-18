#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=1G

start=$(date +%s)

module load fastqc
fastqc -o qc raw_data/*.fastq.gz

end=$(date +%s)
runtime=$((end - start))

echo "Finished"
echo "$runtime"
echo "$(date)"
