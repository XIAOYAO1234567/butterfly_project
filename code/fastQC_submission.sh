######## Quality check ########
#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 1
#$ -l h_rt=240:0:0
#$ -l h_vmem=1G

# Starting time
start=$(date +%s)

# Load FastQC module
module load fastqc

# Run FastQC on all .fastq files in raw_data, output to qc/
fastqc -o qc raw_data/*.fastq.gz

#Finishing timestamp
end=`date +%s`
runtime=$((end - start))

echo Finished
echo $runtime
echo $(date)

