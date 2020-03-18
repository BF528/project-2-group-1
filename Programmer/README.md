This folder has the following contents

qsub scripts to run commands for:
- TopHat
- RseQC utilities (geneBody_coverage, inner_distance, bam_stat) 
- Cufflinks
- cuffdiff

The following outputs 
- bam_stat.txt (has the mapping metrics)
- geneBodymm9.geneBodyCoverage.curves.pdf ( line graph for read coverage)
- innerdistancemm9.inner_distance_plot.pdf (histogram plot of inner distance densities)
- results.txt ( resulats of samtools flagstat for alignment metrics)
- FPKMplot.png ( bar plot of gene FPKM values)
- genes.fpkm_tracking (FPKM values that were given by cufflinks)
