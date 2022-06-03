#############################################################################################
# Figure S7A: Correlation of TuORF normalized read counts vs TuORF transcript mRNA levels
#############################################################################################
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

TuORF_gene <- table(ORF_max_filt_uORF$gene_id)
TuORF_gene_df <- data.frame(gene_id=names(TuORF_gene),freq=as.numeric(TuORF_gene))
head(TuORF_gene_df)
TuORF_gene_df1 <- TuORF_gene_df %>% filter(freq==1)
TuORF_gene_df2 <- TuORF_gene_df %>% filter(freq>1)

#remove genes with more than one uORF from TE_CDS (so others will not include them)
TE_CDS <- TE_CDS %>% filter(!(gene_id %in% TuORF_gene_df2$gene_id))

ORF_max_filt_uORF <- ORF_max_filt_uORF %>% filter(gene_id %in% TuORF_gene_df1$gene_id)

ORF_max_filt_mORF <- ORF_max_filt %>% filter(transcript_id %in% ORF_max_filt_uORF$transcript_id,category=="ORFs_ccds")%>% group_by(transcript_id) %>% top_n(n=1, wt = ORF_length)
ORF_max_filt_mORF <- ORF_max_filt_mORF %>% mutate(mORF_P_sites=ORF_P_sites,mORF_length=ORF_length) %>% dplyr::select(transcript_id,mORF_P_sites,mORF_length)

ORF_max_filt_uORF_mORF <- inner_join(ORF_max_filt_uORF,ORF_max_filt_mORF) %>% mutate(uORF_mORF_Ribo_ratio=(ORF_P_sites/ORF_length)/(mORF_P_sites/mORF_length)) %>% 
  mutate(TuORF_TE=(ORF_P_sites/ORF_RNA_sites)) %>% mutate(TuORF_norm_read_count=ORF_P_sites/ORF_length)

dim(ORF_max_filt_uORF_mORF) #[1] 1438   42

ORF_max_filt_uORF_mORF2 <- left_join(ORF_max_filt_uORF_mORF,TE_CDS,by="gene_id")

F2I <- ggplot(ORF_max_filt_uORF_mORF2, aes(y=log2(TuORF_norm_read_count), x=log2(RNA),alpha=0.001)) + 
  geom_point()+
  geom_smooth(method=lm)+ theme_classic()
F2I
cor(log2(ORF_max_filt_uORF_mORF2$TuORF_norm_read_count), log2(ORF_max_filt_uORF_mORF2$RNA),use="pairwise.complete.obs")
#0.38
cor.test(log2(ORF_max_filt_uORF_mORF2$TuORF_norm_read_count), log2(ORF_max_filt_uORF_mORF2$RNA),use="pairwise.complete.obs")
#p-value < 2.2e-16

#############################################################################################
# Figure S7B: Expression levels for TuORFs, UuORFs and no uORF genes
#############################################################################################
#rm(list=ls())
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

###Find sequence predicted from MAX expressed isoforms
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/Desktop/uORFs_miRNA/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
threeUTRByTx <- threeUTRsByTranscript(txdb,use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 22136
fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,9)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])))
fiveUTR_ORFs_gene_names <- substr(fiveUTR_ORFs_tx_names,1,9)

# LOad RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")

Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id)]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id

RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
TE <- Ribo$TPM/RNA$TPM
TE_CDS <- data.frame(gene_id=RNA$gene_id,
                     RNA=RNA$TPM,
                     Ribo=Ribo$TPM,
                     TE=TE)
# Remove not translated transcripts
TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)
TE_CDS <- TE_CDS %>% filter(TE>0 & !is.infinite(TE))

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq-uORFs",ifelse(TE_CDS$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted-uORFs", "Others"))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("Riboseq-uORFs", "Predicted-uORFs", "Others"), labels=c("Riboseq-uORFs", "Predicted-uORFs", "Others"))
TE_CDS %>% group_by(Category) %>% dplyr::summarise_at(vars(-gene_id),list(mean = mean, median = median))


TE_CDS2 <- TE_CDS %>% filter(Category %in% c("Riboseq-uORFs", "Predicted-uORFs")) %>% droplevels()
TE_CDS2$Category <- factor(TE_CDS2$Category, levels=c("Predicted-uORFs","Riboseq-uORFs"), labels=c("UuORFs","TuORF"))

library(plyr)
mu <- ddply(TE_CDS2, "Category", summarise, grp.mean=median(RNA))

pTE_CDS_hist <- ggplot(TE_CDS2,aes(x=RNA,color=Category,fill=Category)) +
  geom_histogram(alpha=0.8, position = 'identity') +
  xlim(c(0,70)) +
  theme_classic() + geom_vline(data=mu, aes(xintercept=grp.mean, color=Category),linetype = 9,size=1)
pTE_CDS_hist


pTE_CDS_fpoly <- ggplot(TE_CDS2,aes(x=RNA,color=Category)) +
  geom_freqpoly(binwidth = 1) +
  xlim(c(0,100)) +
  theme_classic() + 
  xlab("RNA (TPM)")+
  ylab("Numbers") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Category),linetype = 2,size=0.5)

pTE_CDS_fpoly

# ggsave("~/Desktop/uORFs_miRNA/Figures_May052021/Figure2D_frepoly_TuORFs_vs_UuORFs.pdf",
#        plot = pTE_CDS_fpoly,
#        units = "in",
#        width = 5,
#        height = 4)
