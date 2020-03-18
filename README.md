# Project Description

This repository contains files related to the project "Transcriptional Profile of Mammalian Cardiac Regeneration with mRNA seq" (https://bf528.readthedocs.io/en/latest/content/projects/project_2_rnaseq_1/project_2_rnaseq_1.html)

# Contributors

Divya (Programmer), Garima (Analyst), Marlene (Data Curator), Xudong (Biologist)

# Repository Contents

Contents of the repository are given below:

Analyst folder: contains Analyst results for the project.

Programmer folder: 
Contains qsub files to run commands for:
- TopHat
- RseQC utilities (geneBody_coverage, inner_distance, bam_stat) 
- Cufflinks
- cuffdiff

Contains the following outputs 
- bam_stat.txt (has the mapping metrics)
- geneBodymm9.geneBodyCoverage.curves.pdf ( line graph for read coverage)
- innerdistancemm9.inner_distance_plot.pdf (histogram plot of inner distance densities)
- results.txt ( resulats of samtools flagstat for alignment metrics)
- FPKMplot.png ( bar plot of gene FPKM values)
- genes.fpkm_tracking (FPKM values that were given by cufflinks)


Data Folder : Contains the data that was used for downstream analysis.
