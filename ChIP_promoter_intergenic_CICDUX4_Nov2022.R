library(data.table)
library(dplyr)


setwd("/Volumes/cuyler/ucsf_okimoto_lab/nick_chipseq_project/hg38_alignments/bams/final")



promoters = fread("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER_promoters.bed")
intergenic = fread("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER_intergenic.bed")

CIC_bulk = fread("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_trunc.bed")

CIC_promoters_accessions = CIC_bulk[CIC_bulk$V4 %in% promoters$V1,]
CIC_intergenic_accessions = CIC_bulk[CIC_bulk$V4 %in% intergenic$V1,]

write.table(CIC_promoters_accessions, "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_trunc_promoters.bed",row.names = F, col.names = F, quote = F, sep = "\t")
write.table(CIC_intergenic_accessions, "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_trunc_intergenic.bed",row.names = F, col.names = F, quote = F, sep = "\t")
