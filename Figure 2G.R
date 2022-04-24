```
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

###################################
# dataset 6: 
# Arabidopsis tissue atlas E-MTAB-7978
# rm(list=ls())
dataset6 <- read.xlsx("~/Desktop/uORFs_miRNA/Figures_May052021/E-MTAB-7978-query-results.tpms.xlsx",sheet = 1,startRow =5)
head(dataset6)
library(tidyr)
dataset6 <- dataset6 %>% replace(is.na(.), 0)
colnames(dataset6)[1] <- "gene_id"
dim(dataset6) #[1] 20897     7
dataset6$TPM_MEAN <- rowMeans(dataset6[,3:ncol(dataset6)])
head(dataset6,20) 
dataset6 <- dataset6 %>% filter(gene_id %in% fiveUTRByGeneNames)
dataset6$Category <- ifelse(dataset6$gene_id %in% ORF_max_filt_uORF_id,"atuORFs",ifelse(dataset6$gene_id %in% seq_defined_uORFs, "squORFs", "Others"))
dataset6$Category <- factor(dataset6$Category,levels=c("atuORFs", "squORFs", "Others"), labels=c("atuORFs", "squORFs","no uORF genes"))

table(dataset6$Category)

pRNA <- ggplot(dataset6, aes(x=TPM_MEAN,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,40)+
  labs(title = "dataset6: Col-0 plants under different growth \n conditions from multiple studies") +
  theme_classic() +
  scale_color_manual(values = c("atuORFs" = "#F8766D", "squORFs"="#00BA38", "no uORF genes" = "#619CFF")) 
pRNA

dataset6_mean_by_group <- aggregate(dataset6[,3:(ncol(dataset6)-1)], list(dataset6$Category), FUN=median) %>% as.data.frame()
head(dataset6_mean_by_group)
colnames(dataset6_mean_by_group)[1] <- "Category"

dataset6_mean_by_group$TPM_MEAN <- rowMeans(dataset6_mean_by_group[,3:ncol(dataset6_mean_by_group)])

# dataset6_mean_by_group <- dataset6[,3:ncol(dataset6)] %>% 
#   group_by(Category) %>%
#   summarise_each(funs(mean))

library(reshape2)
dataset6_mean_by_group_melt <- melt(dataset6_mean_by_group)

dp <- ggplot(dataset6_mean_by_group_melt,aes(x=variable,fill=Category,y=value))
dp <- dp + geom_dotplot(binaxis="y")
print(dp)

dataset6m <- melt(dataset6)
p<- ggplot(dataset6m, aes(x=variable, y=value, color=Category)) + ylim(0,85) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab("Samples") + ylab("RNA (TPM)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

#ggsave("~/Desktop/uORFs_miRNA/Figures_May052021/FigureX_Boxplot_Arabidopsis_tissues_RNA-seq_TPM.pdf",
       plot = p,
       units = "in",
       width = 13,
       height = 6)

```
