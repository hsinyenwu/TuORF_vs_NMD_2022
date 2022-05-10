'''
rm(list=ls())
library(openxlsx)
library(dplyr)
library(reshape)
library(ggplot2)
library(ggtext)
library(grid)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(EnvStats)

#Load 5UTR containing gene list
fiveUTRByGeneNames <- read.delim("~/Desktop/uORFs_miRNA/genes_w_5'UTR.tsv",header=T)
fiveUTRByGeneNames <- fiveUTRByGeneNames$gene_id

# Load NMD targets
NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
colnames(NMD) <- "gene_id"

#Load datasets
RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
head(ORF_max_filt_uORF,2)
nrow(ORF_max_filt_uORF)

ORF_max_filt_uORF <- ORF_max_filt_uORF %>% mutate(pept_length=nchar(ORF_pept))

# Only analyze 1 uORF genes
atuORF_gene <- table(ORF_max_filt_uORF$gene_id)
atuORF_gene_df <- data.frame(gene_id=names(atuORF_gene),freq=as.numeric(atuORF_gene))
head(atuORF_gene_df)
atuORF_gene_df1 <- atuORF_gene_df %>% filter(freq==1)
atuORF_gene_df2 <- atuORF_gene_df %>% filter(freq>1)

#remove genes with more than one uORF from TE_CDS (so others will not include them)
#uORF_NC is length normalized count of TuORFs
ORF_max_filt_uORF <- ORF_max_filt_uORF %>% filter(gene_id %in% atuORF_gene_df1$gene_id) %>% mutate(uORF_NC=ORF_P_sites/ORF_length)

ORF_max_filt_uORF_RNA <- inner_join(ORF_max_filt_uORF,RNA,by="gene_id") 
ORF_max_filt_uORF_RNA <- ORF_max_filt_uORF_RNA %>% filter(TPM>0)

library(ggplot2)

# Add the regression line
ggplot(ORF_max_filt_uORF_RNA, aes(x=log2(TPM), y=log2(uORF_NC))) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()  

cor(log2(ORF_max_filt_uORF_RNA$TPM),log2(ORF_max_filt_uORF_RNA$uORF_NC)) #R=0.3994119 ~0.4
'''
