rm(list=ls())
library(dplyr)
library(Biostrings)
library(ORFik)
library(GenomicFeatures)
library(Rsamtools)

FA <- FaFile("/mnt/research/riboplant/Reference/TAIR10_chr_all_2.fas")

ORF_max_filt <- read.delim(file="/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/RiboTaper_CTRL_merged/ORFs_max_filt",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")

txdb_max <- makeTxDbFromGFF("/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/max_exp_isoform_gtf/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByGene_max <- exonsBy(txdb_max,by='gene')
exonByTx_max <- exonsBy(txdb_max,by='tx',use.names=T)
fiveUTRByTx_max <- fiveUTRsByTranscript(txdb_max,use.names=T)
fiveUTR_seqs_max <- extractTranscriptSeqs(FA,fiveUTRByTx_max)

# Calculate MFE for main ORF Kozak using max expressed isoforms
ORF_max_filt2 <- ORF_max_filt %>% 
  filter(transcript_id %in% names(exonByTx_max), !(transcript_id %in% ORF_max_filt_uORF$transcript_id)) %>% 
  filter(annotated_start>60) %>% 
  filter(category=="ORFs_ccds")

nrow(ORF_max_filt2) #15540
#Remove duplicated gene_id to prevent double/more count
ORF_max_filt2 <- ORF_max_filt2[!duplicated(ORF_max_filt2$gene_id),]
nrow(ORF_max_filt2) #15166

#calculate nucleotide frequency (collapse all sequences)
#annotated_start_nt_seq <- c()
#gc()
#for(i in seq_along(ORF_max_filt2$gene_id)){
#  if (i%%1000==0) print(i)
#  tx_sequence <- as.character(unlist(getSeq(FA,exonByTx_max[[ORF_max_filt2$transcript_id[i]]])))
#  aORF_start <- ORF_max_filt2$annotated_start[i]
#  annotated_start_nt_seq[i] <- substr(tx_sequence,aORF_start-60,aORF_start+49)
#}
#
#names(annotated_start_nt_seq) <- ORF_max_filt2$transcript_id
#annotated_start_nt_seq <- annotated_start_nt_seq[nchar(annotated_start_nt_seq)==110]
#head(annotated_start_nt_seq)
#length(annotated_start_nt_seq) #15166
#
#save(annotated_start_nt_seq,file="/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/RNA_MFE/annotated_start_nt_seq.Rda")
load("/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/RNA_MFE/annotated_start_nt_seq.Rda")

#Loop through all rows
MFE_df_annotated_ORF_all <- data.frame(pos=seq(-43,31,1))

# 1:15166
for(i in 1:15166){
  if(i %% 100==0) print(i)
  tx_sequence <- as.character(unlist(getSeq(FA,exonByTx_max[[names(annotated_start_nt_seq)[i]]])))
  line <- annotated_start_nt_seq[i]
  lines <- sapply(1:75,function(x) substr(line,x,x+35)) 
  file.create("uORF_plus_minus50.txt")
  write(lines,file="uORF_plus_minus50.txt")
  file.remove("output")
  system("/mnt/home/larrywu/anaconda3/bin/RNAfold -i uORF_plus_minus50.txt > output")
  readFile <- scan("output", what = character(), sep = "\n", quiet = TRUE)
  MFE <- sapply(seq(1,150,2),function(x) as.numeric(substr(readFile[[x+1]], nchar(readFile[[x]])+4, nchar(readFile[[x+1]])-1)))
  MFE_df <- data.frame(MFE=MFE)
  colnames(MFE_df) <- names(annotated_start_nt_seq)[i]
  MFE_df_annotated_ORF_all <- cbind(MFE_df_annotated_ORF_all,MFE_df)
  file.remove("uORF_plus_minus50.txt")
  file.remove("output")
}

ncol(MFE_df_annotated_ORF_all)
MFE_df_annotated_ORF_all[1:10,1:10]

save(MFE_df_annotated_ORF_all,file="/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/RNA_MFE/MFE_df_annotated_ORF_all_uORF_starts.Rda")

