library(data.table)
library(ChIPseeker)
library(dplyr)
library(stringr)
library(RColorBrewer)


setwd("/Volumes/cuyler/ucsf_okimoto_lab/nick_chipseq_project/hg38_alignments/bams/final")

#load in data

peaks = readPeakFile("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc.bed")

HOMER_peaks = fread(input = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER.txt")

colnames(HOMER_peaks)[1] = "peak_name"


#generate coverage plot for CIC peaks
pdf(file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_covplot.pdf",width = 8, height = 8)
covplot(peaks)
dev.off()


#generate pie chart of annotations for CIC peak genomic locations
HOMER_peaks_truncAnno = mutate(HOMER_peaks, truncAnnotation = str_extract(Annotation, "[[:alpha:]]+"))

HOMER_peaks_truncAnno$truncAnnotation[HOMER_peaks_truncAnno$truncAnnotation=="non" | HOMER_peaks_truncAnno$truncAnnotation=="TTS"] = "Other"

truncAnno_counts=as.numeric(as.matrix(table(HOMER_peaks_truncAnno$truncAnnotation)))

pdf(file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER_truncAnno_pie.pdf", width = 8, height = 8)
names(truncAnno_counts) = c("Exon", "Intergenic", "Intron", "Other", "Promoter", "UTR")
percents = round(100*(truncAnno_counts/sum(truncAnno_counts)),1)
pie(truncAnno_counts, labels=percents, main = "CIC Peak Genomic Locations", col = brewer.pal(length(truncAnno_counts),"Pastel1"), cex = 0.75)
legend("topright", names(truncAnno_counts), cex = 0.8, fill = brewer.pal(length(truncAnno_counts),"Pastel1"))
dev.off()

#export list of unique annotated closest gene symbols for use in gene ontology analysis
write.table(unique(HOMER_peaks$`Gene Name`), file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER_unique_genes.txt",row.names = F, col.names = F, quote = F)
