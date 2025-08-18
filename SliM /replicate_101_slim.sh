#!/bin/bash
#$ -cwd
#$ -l h_vmem=12G
#$ -pe smp 1
#$ -l h_rt=8:00:00
#$ -j y

source ~/miniforge3/bin/activate slim_env

rep=$1
outdir="SLiM_update/outputs/NC_069101.1/replicate/rep${rep}"
mkdir -p "$outdir"

slim -d outdir='"'"$outdir"'"' simulate_NC_101.slim
