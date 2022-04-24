```
###################################
# dataset 7: 
#tomato 
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
```


