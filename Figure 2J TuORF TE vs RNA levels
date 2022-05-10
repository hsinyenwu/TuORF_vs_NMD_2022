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

TE_CDS <- data.frame(gene_id=RNA$gene_id,
                     RNA=RNA$TPM,
                     Ribo=Ribo$TPM,
                     TE=Ribo$TPM/RNA$TPM)
# Remove not translated transcripts
TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)
# TE_CDS <- TE_CDS %>% filter(TE>0 & !is.infinite(TE))
TE_CDS <- TE_CDS %>% filter(RNA>0)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
head(ORF_max_filt_uORF,2)
nrow(ORF_max_filt_uORF)

ORF_max_filt_uORF <- ORF_max_filt_uORF %>% mutate(pept_length=nchar(ORF_pept))

atuORF_gene <- table(ORF_max_filt_uORF$gene_id)
atuORF_gene_df <- data.frame(gene_id=names(atuORF_gene),freq=as.numeric(atuORF_gene))
head(atuORF_gene_df)
atuORF_gene_df1 <- atuORF_gene_df %>% filter(freq==1)
atuORF_gene_df2 <- atuORF_gene_df %>% filter(freq>1)

#remove genes with more than one uORF from TE_CDS (so others will not include them)
TE_CDS <- TE_CDS %>% filter(!(gene_id %in% atuORF_gene_df2$gene_id))

ORF_max_filt_uORF <- ORF_max_filt_uORF %>% filter(gene_id %in% atuORF_gene_df1$gene_id)

#ORF_P_sites<-Total sum of P-sites position mapped to the ORF
#ORF_RNA_sites<-Total sum of RNA-sites position mapped to the ORF
#P_sites_sum<-Total sum of P-sites position mapped to the transcript
#RNA_sites<-Total sum of RNA-sites position mapped to the transcript

ORF_max_filt_mORF <- ORF_max_filt %>% filter(transcript_id %in% ORF_max_filt_uORF$transcript_id,category=="ORFs_ccds")%>% group_by(transcript_id) %>% top_n(n=1, wt = ORF_length)
length(unique(ORF_max_filt_uORF$transcript_id))
length(unique(ORF_max_filt_mORF$transcript_id))
ORF_max_filt_mORF <- ORF_max_filt_mORF %>% mutate(mORF_P_sites=ORF_P_sites,mORF_length=ORF_length) %>% dplyr::select(transcript_id,mORF_P_sites,mORF_length)

ORF_max_filt_uORF_mORF <- inner_join(ORF_max_filt_uORF,ORF_max_filt_mORF) %>% mutate(uORF_mORF_Ribo_ratio=(ORF_P_sites/ORF_length)/(mORF_P_sites/mORF_length)) %>% 
  mutate(TuORF_TE=(ORF_P_sites/ORF_RNA_sites))

dim(ORF_max_filt_uORF_mORF) #[1] 1438   42
ORF_max_filt_uORF_low <- dplyr::top_n(ORF_max_filt_uORF_mORF,-719,uORF_mORF_Ribo_ratio)
ORF_max_filt_uORF_high <- dplyr::top_n(ORF_max_filt_uORF_mORF,719,uORF_mORF_Ribo_ratio)

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% NMD$gene_id,"NMD",ifelse(TE_CDS$gene_id %in% ORF_max_filt_uORF_low$gene_id, "Low", ifelse(TE_CDS$gene_id %in% ORF_max_filt_uORF_high$gene_id,"High","Others")))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("NMD", "Low", "High", "Others"), labels=c("NMD", "Low", "High", "Others"))

TE_CDS %>% group_by(Category) %>% dplyr::dplyr::summarise_at(vars(-gene_id),list(mean = mean, median = median))

ORF_max_filt_uORF_mORF2 <- left_join(ORF_max_filt_uORF_mORF,TE_CDS,by="gene_id")

library(ggplot2)

# Add the regression line
ggplot(ORF_max_filt_uORF_mORF2, aes(y=log2(RNA), x=log2(uORF_mORF_Ribo_ratio),alpha=0.05)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()  

ggplot(ORF_max_filt_uORF_mORF2, aes(y=log2(TE), x=log2(uORF_mORF_Ribo_ratio),alpha=0.05)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()

ggplot(ORF_max_filt_uORF_mORF2, aes(x=log2(TE), y=log2(RNA),alpha=0.05)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()

ggplot(ORF_max_filt_uORF_mORF2, aes(y=log2(RNA), x=log2(TuORF_TE),alpha=0.05)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()  

ggplot(ORF_max_filt_uORF_mORF2, aes(y=log2(TE), x=log2(TuORF_TE),alpha=0.05)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()
