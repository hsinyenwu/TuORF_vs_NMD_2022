```
#########################
rm(list=ls())
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

#Seq-predicted uORFs
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/Desktop/uORFs_miRNA/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 22136

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])))
fiveUTR_ORFs_gene_names <- substr(fiveUTR_ORFs_tx_names,1,9)
length(unique(fiveUTR_ORFs_gene_names)) #[1] 8122

#RiboTaper predicted uORFs
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
ORF_max_filt_uORF_id <- ORF_max_filt_uORF$gene_id
seq_defined_uORFs <- setdiff(fiveUTR_ORFs_gene_names,ORF_max_filt_uORF$gene_id)

length(unique(fiveUTR_ORFs_gene_names))
length(unique(seq_defined_uORFs))
#list of 5'UTR containing genes
fiveUTRByGeneNames <- read.delim("~/Desktop/uORFs_miRNA/genes_w_5'UTR.tsv",header=T)
fiveUTRByGeneNames <- fiveUTRByGeneNames$gene_id

# dataset5:RNA-seq seedling, root, and floral buds of 17 Arabidopsis thaliana accessions
dataset5 <- read.xlsx("~/Desktop/uORFs_miRNA/Figures_May052021/E-GEOD-53197-query-results.tpms.xlsx",sheet = 1,startRow =5)
head(dataset5)
library(tidyr)
dataset5 <- dataset5 %>% replace(is.na(.), 0)
colnames(dataset5)[1] <- "gene_id"
dim(dataset5) #[1] 20897     7
dataset5$TPM_MEAN <- rowMeans(dataset5[,3:ncol(dataset5)])
head(dataset5,20) 
dataset5 <- dataset5 %>% filter(gene_id %in% fiveUTRByGeneNames)
dataset5$Category <- ifelse(dataset5$gene_id %in% ORF_max_filt_uORF_id,"atuORFs",ifelse(dataset5$gene_id %in% seq_defined_uORFs, "squORFs", "Others"))
dataset5$Category <- factor(dataset5$Category,levels=c("atuORFs", "squORFs", "Others"), labels=c("atuORFs", "squORFs","no uORF genes"))

table(dataset5$Category)

dataset5_mean_by_group <- aggregate(dataset5[,3:(ncol(dataset5)-1)], list(dataset5$Category), FUN=median) %>% as.data.frame()
head(dataset5_mean_by_group)
colnames(dataset5_mean_by_group)[1] <- "Category"

dataset5_mean_by_group$TPM_MEAN <- rowMeans(dataset5_mean_by_group[,3:ncol(dataset5_mean_by_group)])

# dataset5_mean_by_group <- dataset5[,3:ncol(dataset5)] %>% 
#   group_by(Category) %>%
#   summarise_each(funs(mean))

library(reshape2)
dataset5_mean_by_group_melt <- melt(dataset5_mean_by_group)

dp <- ggplot(dataset5_mean_by_group_melt,aes(x=variable,fill=Category,y=value))
dp <- dp + geom_dotplot(binaxis="y")
print(dp)

dataset5m <- melt(dataset5)
p<- ggplot(dataset5m, aes(x=variable, y=value, color=Category)) + ylim(0,85) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab("Samples") + ylab("RNA (TPM)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p
```
