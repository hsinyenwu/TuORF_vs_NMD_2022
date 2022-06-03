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
  xlim(c(0,5000)) +
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
