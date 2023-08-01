###
# This script is attempting to look at DNA motifs occurring in CIC-DUX4 peaks
# Cuyler Luck
# Contact: cuyler.luck@ucsf.edu or ross.okimoto@ucsf.edu

#load required libraries
library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(spgs)
library(pheatmap)
library(ggrepel)

#set the working directory
setwd("/Volumes/cuyler/ucsf_okimoto_lab/nick_chipseq_project/hg38_alignments/bams/final")

#read in FASTA sequences of high-confidence CIC-DUX4 peaks (these were generated during earlier ChIP-seq analysis)
#these are already reduced to unique genomic regions only (i.e. identical genomic regions that were present in several significant summits are only represented once)
peaks = fread(file = "CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_trunc_uniqFourth_unique_fullPeak_FASTA.fasta", header = F, quote = F)

#this function will make all possible 8-mers from a sequence, count them in a table, and return that table
getFracRepeatKmer = function(sequence){
  #I want R to make a vector containing every 8-mer (overlapping OK) in this sequence.
  # AKA all subsequences of length 8 beginning with the first one and moving the starting base forward by one until not possible anymore
  # So the end should be the length - 7, and substring is from the counter to counter + 7
  #make the sequence uppercase first, this is not always true due to the nature of the reference FASTA file (soft masked)
  sequence = toupper(sequence)
  
  kmers = c()
  for(counter in seq(1, nchar(sequence)-7, 1)){
    kmers = c(kmers, substring(sequence, counter, counter + 7))
  } 
  
  
  
  #now we can get a table of how frequently each unique kmer occurs
  #and how frequently per total # of kmers in the sequence (divide instance of kmer by length of kmers vector)
  kmer_table = as.data.frame(table(kmers))
  kmer_table = dplyr::mutate(kmer_table, freq.per.kmer = Freq/length(kmers), freq = Freq)
  
  #return the kmer table
  return(kmer_table)
}


#this code cycles through all the sequences in peaks and runs the above function on them, pulling out a few parameters:
#1) the kmer table for that sequence
#2) the name of the peak for that sequence
#3) the actual sequence
kmer_tables = list()
seq_name = c()
sequences = c()

for(num in seq(1,length(peaks$V1),2)){
  output = getFracRepeatKmer(peaks[num+1]$V1) #run getFracRepeatKmer for the sequence corresponding to the peak name in row $num
  kmer_tables = append(kmer_tables, output) #add the kmer_table to the master list of kmer_tables
  norm.freq = data.frame(output$kmers,output$freq.per.kmer) #make a temporary df to hold the kmers and their frequency per kmer
  raw.freq = data.frame(output$kmers,output$freq) #make a temporary df to hold the kmers and their raw frequency
  colnames(norm.freq) = c("kmer", substring(peaks[num]$V1,2)) #change colnames of these temp dfs to give the norm./raw frequencies a name based on sequence name (excluding ">")
  colnames(raw.freq) = c("kmer", substring(peaks[num]$V1,2))
  seq_name = c(seq_name, peaks[num]$V1) #add sequence name to the master list of sequence names
  if(!exists("norm.freq.df")){ #if it doesn't exist already, make a master dataframe for normalized frequencies and for raw frequencies
    norm.freq.df = norm.freq
    raw.freq.df = raw.freq
  }
  else{
    norm.freq.df = merge(norm.freq.df, norm.freq, by = "kmer", all = T) #if it does exist, just merge norm.freq/raw.freq with the master list, keeping all kmers
    raw.freq.df = merge(raw.freq.df, raw.freq, by = "kmer", all = T)
  }
  sequences = c(sequences, toupper(peaks[num+1]$V1)) #add the current sequence to the list of sequences, making uppercase again for good measure
}

#put names and sequences into one dataframe for easier use
sequences_and_names = data.frame(seq_name, sequences)

#need to convert values of NA in norm.freq.df/raw.freq.df to zeroes (this just means that that specific kmer appeared zero times in the given sequence)

norm.freq.df[is.na(norm.freq.df)] = 0
raw.freq.df[is.na(raw.freq.df)] = 0

#write these to files
write.table(norm.freq.df, file = "norm.frac.df.tsv", sep = "\t", quote = F, row.names = F)
write.table(raw.freq.df, file = "raw.freq.df.tsv", sep = "\t", quote = F, row.names = F)


#put kmer as rowname
norm.freq.df.rownames = tibble::column_to_rownames(norm.freq.df,"kmer")
raw.freq.df.rownames = tibble::column_to_rownames(raw.freq.df,"kmer")


#using norm.freq.df.rownames for now (normalized, frequency per kmer of sequence)
#extract kmers that follow the pattern "TGNNTGNN" only (or the reverse complement, "NNCANNCA")

TGNNTGNN_only = norm.freq.df.rownames[str_detect(rownames(norm.freq.df.rownames), "TG..TG..") | str_detect(rownames(norm.freq.df.rownames), "..CA..CA"),]

#calculate row means for each of these TGNNTGNN kmers -- that is, mean frequency per kmer of a given TGNNTGNN kmer across all sequences
rowmeans = rowMeans(TGNNTGNN_only)

#make a dataframe for easier manipulation
rowmeans.df = data.frame(names(rowmeans), rowmeans)

#Now for ease of plotting I want to collapse reverse complements with their cognate kmers
#e.g. TGAATGAA and TTCATTCA form a single value instead of two distinct ones
#to do this I will simply average the rowmeans between the kmers in the kmer:reverse complement pair
#since the number of values comprising the mean is identical for each individual kmer, this is the same thing as if I had calculated the means
#of both kmers together in the first place

already_done = c()
kmer_combined = c()
freq_combined = c()
for(kmer in rowmeans.df$names.rowmeans.){ #for each of the kmers in the names.rowmeans. variable of rowmeans.df
  if(!kmer %in% already_done){ #if that kmer hasn't already been used either on its own or as a reverse complement 
    temp_freq = rowmeans.df[rowmeans.df$names.rowmeans. == kmer,]$rowmeans #get the rowmean for the current kmer
    if(toupper(reverseComplement(kmer)) %in% rowmeans.df$names.rowmeans.){ #check if r.c. is in the list as well. this is necessary because the r.c. is not always in the original kmer list -- will cause issues if dont check
      temp_freq = (temp_freq + rowmeans.df[rowmeans.df$names.rowmeans. == toupper(reverseComplement(kmer)),]$rowmeans)/2 #if the r.c. is in the list, then average its rowmean with the current kmer's rowmean
    } 
    else{ #if the reverse complement wasn't in the list of kmers, then the current kmer should have its current rowmean value cut in half (the r.c. had no instances in any sequence, so rowmean of r.c. = 0)
      temp_freq = temp_freq/2
    }    

    already_done = c(already_done, kmer, toupper(reverseComplement(kmer))) #add the current kmer and its reverse complement to "already_done" so they won't be used again in the loop
    kmer_combined = c(kmer_combined, paste(kmer, " and ", toupper(reverseComplement(kmer)), sep = "")) # put the kmer and reverse complement into a vector
    freq_combined = c(freq_combined, temp_freq) #put the combined frequency into a vector
  } #this loops until all kmers in rowmeans.df$names.rowmeans. have been used, either on their own or as a reverse complement
}

rowmeans.df.combined = data.frame(kmer_combined, freq_combined) #combine kmer_combined and freq_combined for simplicity

#and we can make sure this worked properly by testing a few
(rowmeans.df[rowmeans.df$names.rowmeans. == "AACACTCA",]$rowmeans + rowmeans.df[rowmeans.df$names.rowmeans. == "TGAGTGTT",]$rowmeans)/2
rowmeans.df.combined[rowmeans.df.combined$kmer_combined == "AACACTCA and TGAGTGTT",]$freq_combined
#the values match

(rowmeans.df[rowmeans.df$names.rowmeans. == "CCCAGACA",]$rowmeans + rowmeans.df[rowmeans.df$names.rowmeans. == "TGTCTGGG",]$rowmeans)/2
rowmeans.df.combined[rowmeans.df.combined$kmer_combined == "TGTCTGGG and CCCAGACA",]$freq_combined
#the values match

#what about for testing one whose RC didn't already show up in the list
rowmeans.df[rowmeans.df$names.rowmeans. == "CGCAGCCA",]$rowmeans
#gives a value because this kmer was in the list
rowmeans.df[rowmeans.df$names.rowmeans. == "TGGCTGCG",]$rowmeans
#gives numeric(0) because it wasn't in the list of kmers, i.e. this kmer never showed up in the original sequences
#so the combined mean frequency should just be half of the frequency for the kmer that did appear
rowmeans.df[rowmeans.df$names.rowmeans. == "CGCAGCCA",]$rowmeans / 2
rowmeans.df.combined[rowmeans.df.combined$kmer_combined == "CGCAGCCA and TGGCTGCG",]$freq_combined
#yes, as expected



#make a dotplot of mean frequency per kmer vs. TGNNTGNN 8-mer, ordered by decreasing mean frequency 
pdf("TGNNTGNN.dotplot.pdf", width = 6, height = 6)
ggplot(rowmeans.df.combined, aes(x = reorder(kmer_combined, -freq_combined), y = freq_combined)) + 
  geom_point(size = 3, color = "gray") + 
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 15)) +
  xlab("TGNNTGNN & Reverse Complement 8-mer pairs") +
  ylab("Mean 8-mer Frequency Per 8-mer") +
  geom_label_repel(data = rowmeans.df.combined[rowmeans.df.combined$freq_combined > 2.5e-04,], aes(label = kmer_combined), nudge_x = 150, size = 4) +
  geom_point(data = rowmeans.df.combined[rowmeans.df.combined$freq_combined > 2.5e-04,], color = "#fc0044", size = 3)
dev.off()

#here is a version with no text labels
pdf("TGNNTGNN.dotplot.nolabels.pdf", width = 6, height = 6)
ggplot(rowmeans.df.combined, aes(x = reorder(kmer_combined, -freq_combined), y = freq_combined)) + 
  geom_point(size = 3, color = "gray") + 
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.title = element_text(size = 15)) +
  xlab("TGNNTGNN & Reverse Complement 8-mer pairs") +
  ylab("Mean 8-mer Frequency Per 8-mer") +
  geom_point(data = rowmeans.df.combined[rowmeans.df.combined$freq_combined > 2.5e-04,], color = "#fc0044", size = 3)
dev.off()




#but it's possible that a kmer shows up here either because it's consistently present at lower levels in many sequences, or because
#it's very high in just a few. we can look at this with a heatmap using just the highly-present kmers.
#I'll use the normalized frequencies to choose highly-present kmers, but then I'll look at raw number of appearances of a kmer in a sequence. 
#not normalized by sequence length. this is just to get an idea of how frequently the kmers occur (just once? twice? many times?) in tangible numbers, not proportions

#first, pull out highly-present kmers
chosen = rowmeans.df.combined[rowmeans.df.combined$freq_combined>2.5e-04,]

#next, make it easier to pull out the individual kmers from kmer-rev.comp pairs in each row
chosen = mutate(chosen, kmer1 = substr(kmer_combined, 1, 8))
chosen = mutate(chosen, kmer2 = substr(kmer_combined, 14, 21))

#now go back to the large dataframe with raw frequency of every kmer for all sequences, and just subset to highly-present kmers
chosen_wide = raw.freq.df.rownames[chosen$kmer2,] #this grabs all rows for kmer2, where present
chosen_wide = rbind(chosen_wide, raw.freq.df.rownames[chosen$kmer1,]) #and this adds all rows for kmer1, where present
chosen_wide = as.data.frame(t(chosen_wide))
#now chosen_wide is a dataframe of all sequences (rows) x the 22 chosen kmers (columns) with values being the frequency/kmer of that kmer in that sequence.


#at this point I also want to combine kmer & reverse complements again

chosen_wide_combineKmer = chosen_wide

already_processed = c()
for(kmer in colnames(chosen_wide_combineKmer)){
  
  if(!kmer %in% already_processed){
    rc = toupper(reverseComplement(kmer))
    pair = paste(kmer,"+",rc)
    if(rc %in% colnames(chosen_wide_combineKmer)){
      chosen_wide_combineKmer = dplyr::mutate(chosen_wide_combineKmer, !!pair := chosen_wide_combineKmer[,which(colnames(chosen_wide_combineKmer) == kmer)] + chosen_wide_combineKmer[,which(colnames(chosen_wide_combineKmer) == rc)])
    }
    else{
      chosen_wide_combineKmer = dplyr::mutate(chosen_wide_combineKmer, !!pair := chosen_wide_combineKmer[,which(colnames(chosen_wide_combineKmer) == kmer)])
    }
    already_processed = c(already_processed, kmer, rc)
  }
}

#and we can test that this worked by asking "are any of (kmer + rc == kmer+rc) not true"
any(!(chosen_wide_combineKmer$TGAATGAA + chosen_wide_combineKmer$TTCATTCA == chosen_wide_combineKmer$`TGAATGAA + TTCATTCA`))
any(!(chosen_wide_combineKmer$TGAGTGAG + chosen_wide_combineKmer$CTCACTCA == chosen_wide_combineKmer$`TGAGTGAG + CTCACTCA`))
#should be FALSE

#or by asking if the sum of the kmer and rc columns equals the sum of the (kmer+rc) column
sum(chosen_wide_combineKmer$TGAGTGAG) + sum(chosen_wide_combineKmer$CTCACTCA) == sum(chosen_wide_combineKmer$`CTCACTCA + TGAGTGAG`)
sum(chosen_wide_combineKmer$TGTGTGTG) + sum(chosen_wide_combineKmer$CACACACA) == sum(chosen_wide_combineKmer$`CACACACA + TGTGTGTG`)
#should be TRUE

#ok now we can remove columns that are uncombined kmers
#doing this by keeping columns that have a + in them
chosen_wide_combineKmer = dplyr::select(chosen_wide_combineKmer, grep("\\+",colnames(chosen_wide_combineKmer)))

#manually reordering the kmer columns to match descending order of the dotplot above
col_order = c("TGGATGGA + TCCATCCA","TGAATGAA + TTCATTCA","TCCATTCA + TGAATGGA","TGGATGAA + TTCATCCA","CTCATTCA + TGAATGAG","TTCACTCA + TGAGTGAA","CACACACA + TGTGTGTG","CTCACTCA + TGAGTGAG","CCCATTCA + TGAATGGG","CACATTCA + TGAATGTG","CTCACACA + TGTGTGAG")
chosen_wide_combineKmer = chosen_wide_combineKmer[,col_order]

#and I would also like to order the rows such that the sequences with the highest total raw # of kmer occurrences are at the top, decreasing to the bottom

raw_sums = chosen_wide_combineKmer
raw_sums = dplyr::mutate(raw_sums, rowsum = rowSums(chosen_wide_combineKmer))
raw_sums = tibble::rownames_to_column(raw_sums)
raw_sums = raw_sums[order(raw_sums$rowsum, decreasing = T),]

chosen_wide_combineKmer = chosen_wide_combineKmer[raw_sums$rowname,]


pdf("TGNNTGNN_rawcounts_heatmap.pdf", width = 6, height = 4)
pheatmap(chosen_wide_combineKmer, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,10,length.out=100)), color = colorRampPalette(c("gray90", "#fc0044"))(100))
#please note that the above command restricts the color scale so any values greater than 10 show up as the color for 10
#without this, smaller numbers are easily washed out
dev.off()



#now I would also like to be able to see the genomic features these sequences are annotated for.
#to do this, I need to connect the sequences to the data called by HOMER

#first, read in HOMER data
homer = fread("CIC_IgG_q0.001_DAC_C3_unique_nonblacklisted_chromosomal_chrLabels_trunc_HOMER.txt")

#this is good but unfortunately the naming scheme for chromosomes that is in my data used accessions, not chrXYZ
#so let's reuse some code I wrote elsewhere to convert the HOMER chromosome names to accessions

#manually pair chrXYZ to accessions used in the reference genome
chrs = c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX")
accessions = c("NC_000001.11","NC_000002.12","NC_000003.12","NC_000004.12","NC_000005.10","NC_000006.12","NC_000007.14","NC_000008.11","NC_000009.12","NC_000010.11","NC_000011.10","NC_000012.12","NC_000013.11","NC_000014.9","NC_000015.10","NC_000016.10","NC_000017.11","NC_000018.10","NC_000019.10","NC_000020.11","NC_000021.9","NC_000022.11","NC_000023.11")
translate_chrs = data.frame(chrs = chrs, accessions = accessions)

#make a temporary variable in homer called "temp" and then change it to be the proper accession for that row
homer = dplyr::mutate(homer, accession = "temp")

for(chr in unique(homer$Chr)){
  homer[homer$Chr == chr,]$accession = translate_chrs[translate_chrs$chrs == chr,]$accessions
}

#now we can make the same format of sequence name in "homer" that is used in the rest of the data
#that is, accession+(start-1)+end
#start-1 is needed I think because of bedtools vs. homer formatting. the combined names don't come out identical otherwise
homer = dplyr::mutate(homer, name = paste(accession,Start-1,End,sep=""))

#and let's reduce this to just unique genomic locations (it currently contains all summits, which are distinct peaks but not necessarily separate genomic locations)
#this should be OK because summits are likely annotated for the same things as each other
#do this by ordering by name, and then keeping non-duplicated only
homer = homer[order(homer$name)]
homer = homer[!duplicated(homer$name),]
#this correctly narrows things down to 374 rows

#I also want to simplify the annotations, they are a little too detailed. I wrote code to do this elsewhere:
homer = mutate(homer, Element = str_extract(Annotation, "[[:alpha:]]+"))
homer$Element[homer$Element=="non" | homer$Element=="TTS"] = "other"

#now let's move Element to its own df, put the name as rownames, and use this for annotation on a heatmap
homer_truncAnno = dplyr::select(homer, name, Element)
homer_truncAnno = tibble::column_to_rownames(homer_truncAnno, "name")

#manually specifying some distinct colors for the annotations
anno_colors = list(Element = c(intron = "#7e0e8c", Intergenic = "#ffbb00", promoter = "#006fc3", other = "#ff8340", exon = "#149c3d", UTR = "#888a8c"))


pdf("TGNNTGNN_rawcounts_heatmap_peakElementAnnotated.pdf", width = 6, height = 4)
pheatmap(chosen_wide_combineKmer, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,10,length.out=100)), color = colorRampPalette(c("gray90", "#fc0044"))(100), annotation_row = homer_truncAnno, annotation_colors = anno_colors)
#please note that the above command restricts the color scale so any values greater than 10 show up as the color for 10
#without this, smaller numbers are easily washed out
dev.off()




#let's look at some examples of individual genes

### ETV4
ETV4_homer = homer[homer$`Gene Name` == "ETV4",]
ETV4 = chosen_wide_combineKmer[ETV4_homer$name,]

pdf("ETV4.pdf", width = 6, height = 2.5)
pheatmap(ETV4, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(ETV4),length.out=100)), color = colorRampPalette(c("gray90", "#fc0044"))(100), display_numbers = T, fontsize_number = 8)
dev.off()


### FOXN3

FOXN3_homer = homer[homer$`Gene Name` == "FOXN3",]
FOXN3 = chosen_wide_combineKmer[FOXN3_homer$name,]

pdf("FOXN3.pdf", width = 6, height = 2.5)
pheatmap(FOXN3, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(FOXN3),length.out=100)), color = colorRampPalette(c("gray90", "#fc0044"))(100), display_numbers = T, fontsize_number = 8)
dev.off()


### SPON2

SPON2_homer = homer[homer$`Gene Name` == "SPON2",]
SPON2 = chosen_wide_combineKmer[SPON2_homer$name,]

pdf("SPON2.pdf", width = 6, height = 2.5)
pheatmap(SPON2, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(SPON2),length.out=100)), color = colorRampPalette(c("gray90", "#fc0044"))(100), display_numbers = T, fontsize_number = 8)
dev.off()


### CDH4
CDH4_homer = homer[homer$`Gene Name` == "CDH4",]
CDH4 = chosen_wide_combineKmer[CDH4_homer$name,]

pdf("CDH4.pdf", width = 6, height = 2.5)
pheatmap(CDH4, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(CDH4),length.out=100)), color = colorRampPalette(c("gray90", "#006fc3"))(100), display_numbers = T, fontsize_number = 8)
dev.off()


### POLE
POLE_homer = homer[homer$`Gene Name` == "POLE",]
POLE = chosen_wide_combineKmer[POLE_homer$name,]

pdf("POLE.pdf", width = 6, height = 2.5)
pheatmap(POLE, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(POLE),length.out=100)), color = colorRampPalette(c("gray90", "#006fc3"))(100), display_numbers = T, fontsize_number = 8)
dev.off()

### SPRED2
SPRED2_homer = homer[homer$`Gene Name` == "SPRED2",]
SPRED2 = chosen_wide_combineKmer[SPRED2_homer$name,]

pdf("SPRED2.pdf", width = 6, height = 2.5)
pheatmap(SPRED2, fontsize=5, cluster_cols = F, cluster_rows = F, show_rownames = F, show_colnames=T, breaks = c(seq(0,max(SPRED2),length.out=100)), color = colorRampPalette(c("gray90", "#006fc3"))(100), display_numbers = T, fontsize_number = 8)
dev.off()

# #fc0044
# #ffbb00
# #006fc3
# #fc558f
# #a2ad57
# #e3ac8d
# #85adff
