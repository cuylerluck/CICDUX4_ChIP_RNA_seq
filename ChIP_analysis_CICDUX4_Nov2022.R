library(data.table)
library(dplyr)

setwd("/Volumes/cuyler/ucsf_okimoto_lab/nick_chipseq_project/hg38_alignments/bams/final")

CIC_original = fread("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal.bed")
CIC_chrs = CIC_original


chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
accessions = c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12","NC_000005.10","NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12","NC_000010.11","NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10","NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11","NC_000021.9","NC_000022.11","NC_000023.11")

translate_chrs = data.frame(chrs = chrs, accessions = accessions)

for (accession in translate_chrs$accessions) {
  
  CIC_chrs[CIC_chrs$V1 == accession,]$V1 = translate_chrs[translate_chrs$accessions == accession,]$chrs
  
  
}

write.table(CIC_chrs, file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels.bed", quote = F, row.names = F, col.names = F, sep = "\t")


CIC_chrs_trunc = select(CIC_chrs, c(1:6))

write.table(CIC_chrs_trunc, file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc.bed", quote = F, row.names = F, col.names = F, sep = "\t")


CIC_original_trunc = select(CIC_original, c(1:6))

write.table(CIC_original_trunc, file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_trunc.bed", quote = F, row.names = F, col.names = F, sep = "\t")








