#############################################################################################
# Figure S8: Freq Poly plot for NMD targets with or without TuORFs
#############################################################################################
rm(list=ls())
library(reshape)
library(dbplyr)
library(reshape)
library(ggplot2)
library(ggtext)
library(ORFik)
library(ggpubr)
library(rstatix)
library(gridExtra)
library(grid)
library(EnvStats)
library(openxlsx)
library(GenomicFeatures)
library(GenomicRanges)

#################################################
# Find predicted uORFs, only from max isoforms
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/Desktop/uORFs_miRNA/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
# txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
# txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206_expressed.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
cdsByTx <- cdsBy(txdb,by='tx',use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 22059

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])))
fiveUTR_ORFs_gene_names <- unique(substr(fiveUTR_ORFs_tx_names,1,9)) #[1] 10112

#####################
# Load RiboTaper output
Exp <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(Exp$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
Exp_uORFs <- Exp %>% filter(category=="uORF")

###########################
# Load NMD targets
NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
colnames(NMD) <- "gene_id"

NMD_target_wo_uORFs <- setdiff(NMD$gene_id,Exp_uORFs$gene_id)
length(NMD_target_wo_uORFs) #[1] 284
NMD_target_w_uORFs <- intersect(NMD$gene_id,Exp_uORFs$gene_id)
length(NMD_target_w_uORFs) #[1] 49
uORFs_not_NMD_targets <- setdiff(Exp_uORFs$gene_id,NMD$gene_id)
length(uORFs_not_NMD_targets) #[1] 1719
Others <- setdiff(Exp$gene_id,c(NMD$gene_id,Exp_uORFs$gene_id))
length(Others) #19461

Exp_uORFs %>% filter(gene_id=="AT1G05160")


Exp_uORFs_NMD <- Exp_uORFs %>% filter(gene_id %in% NMD_target_w_uORFs) %>% mutate(ORF_pept_length =nchar(ORF_pept),Category="NMD") %>% dplyr::select(gene_id,ORF_pept_length,Category)
Exp_uORFs_not_NMD <- Exp_uORFs %>% filter(!(gene_id %in% NMD_target_w_uORFs))%>% mutate(ORF_pept_length =nchar(ORF_pept),Category="non_NMD") %>% dplyr::select(gene_id,ORF_pept_length,Category)

Exp_uORFs_if_NMD <- rbind(Exp_uORFs_NMD,Exp_uORFs_not_NMD)

Exp_uORFs_if_NMD$Category <- factor(Exp_uORFs_if_NMD$Category, levels=c("NMD", "non_NMD"), labels=c("NMD", "non_NMD"))

bxp <- ggboxplot(Exp_uORFs_if_NMD, 
                 x = "Category", y = "ORF_pept_length", 
                 color = "Category")+
  ylab("ORF_pept_length")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

bxp

library(plyr)
mu <- ddply(Exp_uORFs_if_NMD, "Category", summarise, grp.mean=median(ORF_pept_length))

fpoly <- ggplot(Exp_uORFs_if_NMD,aes(x=ORF_pept_length,color=Category)) +
  geom_freqpoly(binwidth = 1) +
  xlim(c(0,200)) +
  theme_classic() + 
  xlab("Peptide length")+
  ylab("Numbers") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Category),linetype = 2,size=0.5)

fpoly

