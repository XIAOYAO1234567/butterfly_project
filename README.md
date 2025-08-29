# butterfly_project
# Detecting Signatures of Natural Selection in Butterfly Genomes (anynana)

This project aims to detect signals of natural selection—including both positive and balancing selection—in the genomes of the butterfly species _Bicyclus anynana_.    
By leveraging population genetic statistics and high-throughput sequencing data, we identify genomic regions potentially under selection.

## Project Goals
Analyze whole-genome resequencing data from 5 _Bicyclus anynana_ individuals.
Compute population-genetic statistics (pi, Tajima’s D, LD) and identify candidate windows.
Train and apply a CNN to classify windows as neutral / sweep / balancing.
Functionally annotate genes overlapping consensus candidate regions.


## Tools and Technologies

This project relies on high-performance computing (HPC) and several bioinformatics tools, including:
Core bioinformatics: BWA(sequence alignment), samtools,bcftools(variant calling and filtering), ANGSD(SFS and theta statistics), BedTools(coverage analysis)
Simulations: SLiM
R (visualization and statistical analysis): tidyverse, data.table, ggplot2
Python (MLand plotting): numpy, pandas, scikit-learn, matplotlib, torch 
HPC: SGE (qsub)；
