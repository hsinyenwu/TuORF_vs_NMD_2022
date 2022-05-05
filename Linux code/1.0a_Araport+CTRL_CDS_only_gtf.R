module purge; module load GCC/9.3.0  OpenMPI/4.0.3 R/4.0.3

mkdir -p /mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/CDS_only_gtf

R

library(dplyr)
rm(list=ls())
Araport_gtf <- read.delim("/mnt/home/larrywu/ABA/20181206/reference/Araport11_20181206.gtf",quote="",header=F,sep="\t",stringsAsFactors = F,skip=0)
#/mnt/home/larrywu/ABA/20181206/reference/Araport11_20181206.gtf
table(Araport_gtf$V3)
Araport_gtf_c <- Araport_gtf[which(Araport_gtf$V3=="CDS"),]
Araport_gtf_c$V3 <- "exon"
Araport_gtf_c$gene_id <- sapply(1:nrow(Araport_gtf_c), function(x) gsub("\"","", gsub("gene_id ", replacement = "",unlist(strsplit(Araport_gtf_c$V9[x],"; "))[1])))
Araport_gtf_c$tx_id <- sapply(1:nrow(Araport_gtf_c), function(x) gsub("\"","", gsub("transcript_id ", replacement = "",unlist(strsplit(Araport_gtf_c$V9[x],"; "))[2])))
head(Araport_gtf_c)
MS_gtf_list <- split(Araport_gtf_c,Araport_gtf_c$gene_id)

gene_rows <- function(x) {
  gene <- x[1,]
  gene$V3 <-  "gene"
  gene$V4 <- min(x$V4)
  gene$V5 <- max(x$V5)
  gene$V9 <- paste0("gene_id \"",x$gene_id[1],"\";")
  tx <- gene
  tx$V3 <- "transcript"
  tx$V9 <- paste0(tx$V9," transcript_id \"",x$tx_id[1],"\";")
  rbind(gene,tx,x)[,1:9]
}

MS_gtf2 <- lapply(MS_gtf_list,function(x) gene_rows(x))
MS_gtf2_df <- do.call(rbind,MS_gtf2)
rownames(MS_gtf2_df) <- NULL
MS_gtf2_df$V8 <- "."
head(MS_gtf2_df,20)
write.table(MS_gtf2_df,"/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/CDS_only_gtf/Araport11_20181206_CDS_only.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("/mnt/home/larrywu/Kyle_CTRL/analysis_Araport11_1st+2nd_Riboseq/CDS_only_gtf/Araport11_20181206_CDS_only.gtf",format="gtf", dataSource="Araport",organism="Arabidopsis thaliana")



#Araport_gtf_ncRNA <- read.delim("/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/CTRL_sORFs_expressed.gtf",header=F,sep="\t",quote="",stringsAsFactors = F,skip=0)
#table(Araport_gtf_ncRNA$V3)
#Araport_gtf_c <- Araport_gtf_ncRNA[which(Araport_gtf_ncRNA$V3=="CDS"),]
#Araport_gtf_c$V3 <- "exon"
#Araport_gtf_c$gene_id <- sapply(1:nrow(Araport_gtf_c), function(x) gsub("\"","", gsub("gene_id ", replacement = "",unlist(strsplit(Araport_gtf_c$V9[x],"; "))[1])))
#Araport_gtf_c$tx_id <- sapply(1:nrow(Araport_gtf_c), function(x) gsub("\"","", gsub("transcript_id ", replacement = "",unlist(strsplit(Araport_gtf_c$V9[x],"; "))[2])))
#head(Araport_gtf_c)
#MS_gtf_list <- split(Araport_gtf_c,Araport_gtf_c$gene_id)
#
#gene_rows <- function(x) {
#  gene <- x[1,]
#  gene$V3 <-  "gene"
#  gene$V4 <- min(x$V4)
#  gene$V5 <- max(x$V5)
#  gene$V9 <- paste0("gene_id \"",x$gene_id[1],"\";")
#  tx <- gene
#  tx$V3 <- "transcript"
#  tx$V9 <- paste0(tx$V9," transcript_id \"",x$tx_id[1],"\";")
#  rbind(gene,tx,x)[,1:9]
#}
#MS_gtf2 <- lapply(MS_gtf_list,function(x) gene_rows(x))
#MS_gtf2_df_nc <- do.call(rbind,MS_gtf2)
#rownames(MS_gtf2_df_nc) <- NULL
#MS_gtf2_df_nc$V8 <- "."
#head(MS_gtf2_df_nc,20)
#MS_gtf2_df_both <- rbind(MS_gtf2_df,MS_gtf2_df_nc)
#
#write.table(MS_gtf2_df_both,"/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/Araport11+CTRL_20181206_expressed_CDS_only.gtf",quote = F, row.names = FALSE,col.names = FALSE, sep="\t")
#
#library(GenomicFeatures)
#txdb <- makeTxDbFromGFF("/mnt/home/larrywu/CTRL_arabidopsis/data/assembledGTF/Araport11+CTRL_20181206_expressed_CDS_only.gtf",format="gtf", dataSource="Araport",organism="Arabidopsis thaliana")

##

