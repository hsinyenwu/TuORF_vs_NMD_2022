#############################################################################################
# Figure S5: mRNA levels of TuORFs, UuORFs, and other genes
#############################################################################################

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

#####################
#S5B
#tomato 
#Transcription profiling by high throughput sequencing of tomato root, leaf, flower (2 stages) and fruit (6 stages).

rm(list=ls())
dataset8 <- read.delim("~/Desktop/uORFs_miRNA/E-MTAB-4813-query-results.tpms.tsv",header = T,skip=4)
dataset8 <- dataset8[270:nrow(dataset8),]
head(dataset8)
library(tidyr)
dataset8 <- dataset8 %>% replace(is.na(.), 0)
colnames(dataset8)[1] <- "gene_id"
dim(dataset8) #[1] 23459    20
dataset8$TPM_MEAN <- rowMeans(dataset8[,3:ncol(dataset8)])
head(dataset8,20) 
# ~/Desktop/Tomato_Next/ORFs_max_filt_Root.tsv
# ~/Desktop/Tomato_Next/ORFs_max_filt_Shoot.tsv
Root_ORFs_max_filt <- read.delim("~/Desktop/Tomato_Next/ORFs_max_filt_Root.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
Root_ORFs_max_filt_uORFs <- Root_ORFs_max_filt %>% filter(category=="uORF")
#genes with 5'UTR
FA <- FaFile("~/Desktop/Tomato/Reference/seqs_srcdir/S_lycopersicum_3.00.fa")
txdb <- makeTxDbFromGFF("~/Desktop/Tomato_paper_version//Tomato_Root_20181118/SupTables/Tomato_Root_ITAG_vs_de_novo_10312018.ixyous+ITAG3.2_rm_start_stop.gtf",format="gtf", dataSource="ITAG3.2",organism="Solanum lycopersicum")

exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 17868

fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,16)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)
dataset8 <- dataset8 %>% filter(gene_id %in% fiveUTRByGeneNames)

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=F,groupByTx=F, minimumLength=9)
length(unique(fiveUTR_ORFs)) #42328
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])[1]))
length(fiveUTR_ORFs_tx_names) #42328
length(unique(fiveUTR_ORFs_tx_names))#9166
fiveUTR_ORFs_gene_names <- substr(unique(fiveUTR_ORFs_tx_names),1,16)
length(fiveUTR_ORFs_gene_names) #9166

tail(fiveUTR_ORFs_gene_names)
sum(grepl("Solyc",fiveUTR_ORFs_gene_names)) #9166
Only_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names %in% Root_ORFs_max_filt_uORFs$gene_id)]

data8 <- dataset8 %>% filter(gene_id %in% fiveUTRByGeneNames)
data8$Category <- ifelse(data8$gene_id %in% Root_ORFs_max_filt_uORFs$gene_id,"Riboseq",ifelse(data8$gene_id %in% Only_predicted_uORFs, "Predicted", "Others"))
table(data8$Category)
# Others Predicted   Riboseq 
# 7123      6422      1088
data8$Category <- factor(data8$Category, levels=c("Riboseq", "Predicted","Others"), labels=c("ATuORFs", "SPuORFs", "Others"))

table(data8$Category)

pRNA <- ggplot(data8, aes(x=TPM_MEAN,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,40)+
  labs(title = "Expression of difernt tissues in tomato") +
  theme_classic() +
  scale_color_manual(values = c("ATuORFs" = "#F8766D", "SPuORFs"="#00BA38", "Others" = "#619CFF")) 
pRNA

data8_mean_by_group <- aggregate(data8[,3:(ncol(data8)-1)], list(data8$Category), FUN=median) %>% as.data.frame()
head(data8_mean_by_group)
colnames(data8_mean_by_group)[1] <- "Category"

data8_mean_by_group$TPM_MEAN <- rowMeans(data8_mean_by_group[,3:ncol(data8_mean_by_group)])

# data8_mean_by_group <- data8[,3:ncol(data8)] %>% 
#   group_by(Category) %>%
#   summarise_each(funs(mean))

library(reshape2)
data8_mean_by_group_melt <- melt(data8_mean_by_group)

dp <- ggplot(data8_mean_by_group_melt,aes(x=variable,fill=Category,y=value))
dp <- dp + geom_dotplot(binaxis="y")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
print(dp)

data8m <- melt(data8)
p<- ggplot(data8m, aes(x=variable, y=value, color=Category)) + ylim(0,85) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab("Samples") + ylab("RNA (TPM)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

ggsave("~/Desktop/uORFs_miRNA/Figures_Aug252021/FigureX_Boxplot_tomato_tissues_RNA-seq_TPM_E-MTAB-4813.pdf",
       plot = p,
       units = "in",
       width = 13,
       height = 6)
############################
# S5C
#tomato dataset 7: 
#Transcription profiling by high throughput sequencing of tomato root, leaf, flower (2 stages) and fruit (6 stages).
"~/Desktop/uORFs_miRNA/Figures_May052021/E-MTAB-4812-query-results.tpms.tsv"# 
rm(list=ls())
#Data download from https://www.ebi.ac.uk/gxa/experiments?kingdom=Plants

dataset8 <- read.delim("~/Desktop/uORFs_miRNA/Figures_May052021/E-MTAB-4812-query-results.tpms.tsv",skip =4)
dataset8 <- dataset8[407:nrow(dataset8),]
head(dataset8)
library(tidyr)
dataset8 <- dataset8 %>% replace(is.na(.), 0)
colnames(dataset8)[1] <- "gene_id"
dim(dataset8) #[1] 27062     7
dataset8$TPM_MEAN <- rowMeans(dataset8[,3:ncol(dataset8)])
head(dataset8,20) 

Root_ORFs_max_filt <- read.delim("~/Desktop/Tomato_Next/ORFs_max_filt_Root.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
Root_ORFs_max_filt_uORFs <- Root_ORFs_max_filt %>% filter(category=="uORF")
#genes with 5'UTR
FA <- FaFile("~/Desktop/Tomato/Reference/seqs_srcdir/S_lycopersicum_3.00.fa")
txdb <- makeTxDbFromGFF("~/Desktop/Tomato_paper_version//Tomato_Root_20181118/SupTables/Tomato_Root_ITAG_vs_de_novo_10312018.ixyous+ITAG3.2_rm_start_stop.gtf",format="gtf", dataSource="ITAG3.2",organism="Solanum lycopersicum")

exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 17868

fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,16)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)
dataset8 <- dataset8 %>% filter(gene_id %in% fiveUTRByGeneNames)

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=F,groupByTx=F, minimumLength=9)
length(unique(fiveUTR_ORFs)) #42328
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])[1]))
length(fiveUTR_ORFs_tx_names) #42328
length(unique(fiveUTR_ORFs_tx_names))#[1] 9166
fiveUTR_ORFs_gene_names <- substr(unique(fiveUTR_ORFs_tx_names),1,16)
length(fiveUTR_ORFs_gene_names) #9166

tail(fiveUTR_ORFs_gene_names)
sum(grepl("Solyc",fiveUTR_ORFs_gene_names)) #9166
Only_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names %in% Root_ORFs_max_filt_uORFs$gene_id)]

data7 <- dataset8 %>% filter(gene_id %in% fiveUTRByGeneNames)
data7$Category <- ifelse(data7$gene_id %in% Root_ORFs_max_filt_uORFs$gene_id,"Riboseq",ifelse(data7$gene_id %in% Only_predicted_uORFs, "Predicted", "Others"))
table(data7$Category)
# Others Predicted   Riboseq 
# 7637      7008      1101 
data7$Category <- factor(data7$Category, levels=c("Riboseq", "Predicted","Others"), labels=c("ATuORFs", "SPuORFs", "Others"))

table(data7$Category)

data7_mean_by_group <- aggregate(data7[,3:(ncol(data7)-1)], list(data7$Category), FUN=median) %>% as.data.frame()
head(data7_mean_by_group)
colnames(data7_mean_by_group)[1] <- "Category"

data7_mean_by_group$TPM_MEAN <- rowMeans(data7_mean_by_group[,3:ncol(data7_mean_by_group)])

library(reshape2)
data7_mean_by_group_melt <- melt(data7_mean_by_group)

data7m <- melt(data7)
p<- ggplot(data7m, aes(x=variable, y=value, color=Category)) + ylim(0,85) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab("Samples") + ylab("RNA (TPM)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

ggsave("~/Desktop/uORFs_miRNA/Figures_Aug252021/FigureX_Boxplot_tomato_tissues_RNA-seq_TPM_E-MTAB-4812.pdf",
       plot = p,
       units = "in",
       width = 13,
       height = 6)

