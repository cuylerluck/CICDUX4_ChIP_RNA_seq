# CICDUX4_ChIP_RNA_seq
Analysis of CIC-DUX4 ChIP-seq and RNA-seq data

This repository houses the code used to analyze ChIP-seq data and RNA-seq data described in Thomas & Luck et al (2023).
These scripts are not plug-and-play, but are meant to provide transparency into how we analyzed the data.
With questions, contact Cuyler Luck at Cuyler.Luck@ucsf.edu or Ross Okimoto at Ross.Okimoto@ucsf.edu.

The files contain code as follows:

# ChIP-seq...

CICDUX4_ChIPseq_documentation.txt: Describes all analysis done for ChIP-seq experiments done with X1 cells that was not done in R. This includes processing/aligning of FASTQ files, peak calling, and extracting sequences for motif analyses.

ChIP_analysis_CICDUX4_Nov2022.R: A simple R script used to interchange chromosome naming between "chrXYZ" and "NC_00000X.YZ". Referenced in the txt documentation.

ChIP_gene_ontology_CICDUX4_Nov2022.R: R script used to make a coverage plot for peaks, generate the piechart of annotated genomic locations, and output a list of gene symbols for gene ontology analysis if desired. Referenced in the txt documentation.

ChIP_promoter_intergenic_CICDUX4_Nov2022.R: Simple R script used to subset peaks to those annotated for promoter and intergenic regions. Referenced in the txt documentation.

motif_analysis.R: Lengthy R script that analyzed 8-mers appearing in high-confidence CIC-DUX4 peaks. A longer description is available in the manuscript, or commented in the code.

# RNA-seq...

documentation_RNAseq_siCIC_X1_CDS2.txt: Describes all non-R analysis done for RNA-seq data of X1/CDS2 cell lines treated with a control or CIC-targeted siRNA. This includes sequence alignment and processing as well as grep commands.

X1_CDS2_siCIC_DE.R: Describes all R analysis done with RNA-seq data from above. This includes edgeR differential expression analyses and lots of plotting.
