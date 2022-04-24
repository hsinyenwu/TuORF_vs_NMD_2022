```
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
```
