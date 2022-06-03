#############################################################################################
# Figure S3 TuORFs vs UuORFs vs Other genes cumulative curve for mRNA levels, TE, mRNA half-lives and protein abundance
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

###Find sequence predicted from MAX expressed isoforms
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
txdb <- makeTxDbFromGFF("~/Desktop/uORFs_miRNA/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 22136

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=F,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]][1])))
length(fiveUTR_ORFs_tx_names) #[1] 17125
fiveUTR_ORFs_gene_names <- substr(fiveUTR_ORFs_tx_names,1,9)

RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
TE <- Ribo$TPM/RNA$TPM

CPuORFs <- read.xlsx("~/Desktop/uORFs_miRNA/TableS2_uORFs_on_CPuORF_containg_5UTR_comprehensive_sing_CPuORF_genes_only_colored.xlsx")

TE_CDS <- data.frame(gene_id=RNA$gene_id,
                     RNA=RNA$TPM,
                     Ribo=Ribo$TPM,
                     TE=TE)
#genes with 5'UTR
fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,9)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)
# write.table(fiveUTRByGeneNames_df,"~/Desktop/uORFs_miRNA/genes_w_5'UTR.tsv",quote = F, row.names = FALSE,col.names = TRUE, sep="\t")
TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)

# Remove not translated transcripts
TE_CDS <- TE_CDS %>% filter(TE>0)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% CPuORFs$gene_id, "CPuORFs",ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(TE_CDS$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others")))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("CPuORFs","Riboseq", "Predicted", "Others"), labels=c("CPuORFs","Riboseq", "Predicted","Others"))

TE_CDS %>% group_by(Category) %>% summarise(median=median(RNA)/7.82)
TE_CDS %>% group_by(Category) %>% summarise(median=median(TE)/1.19)

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(TE_CDS$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others"))
table(TE_CDS$Category)
# CPuORFs    Others Predicted   Riboseq 
# 79     13060      6330      1683 
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))


# mRNA half-lives
half_lives <- read.delim(file="/Volumes/GoogleDrive/My Drive/CPuORF_paper/data-for-analysis/RNA-half-lives-5EU.txt",header=T,stringsAsFactors=F,sep="\t")
half_lives <- half_lives %>% dplyr::select("Gene.ID","Mean")
colnames(half_lives) <- c("gene_id","Mean")
half_lives <- half_lives %>% filter(gene_id %in% fiveUTRByGeneNames)
# half_lives$Category <- ifelse(half_lives$gene_id %in% CPuORFs$gene_id, "CPuORFs",ifelse(half_lives$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(half_lives$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others")))
half_lives$Category <- ifelse(half_lives$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(half_lives$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others"))
half_lives$Category <- factor(half_lives$Category, levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))

# Quantitative peoteomics for Shoot
shoot_proteomics <- read.xlsx("~/Desktop/uORFs_miRNA/Justin_Proteomics2017/data_excel/pmic12946-sup-0004-tables7.xlsx",sheet = 3)
head(shoot_proteomics$`MS/MS.count`)
shoot_proteomics$gene_id <- substr(shoot_proteomics$Protein.IDs,1,9)
shoot_proteomics2 <- shoot_proteomics %>% dplyr::select(gene_id,`MS/MS.count`,`Sequence.length`)
colnames(shoot_proteomics2) <-c("gene_id","MS_count","Pept_length")
shoot_proteomics2 <- shoot_proteomics2 %>% filter(gene_id %in% fiveUTRByGeneNames)
shoot_proteomics2$norm_MS_count <- shoot_proteomics2$MS_count/shoot_proteomics2$Pept_length
# shoot_proteomics2$Category <- ifelse(shoot_proteomics2$gene_id %in% CPuORFs$gene_id, "CPuORFs",ifelse(shoot_proteomics2$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(shoot_proteomics2$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others")))
shoot_proteomics2$Category <- ifelse(shoot_proteomics2$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(shoot_proteomics2$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others"))
shoot_proteomics2$Category <- factor(shoot_proteomics2$Category,levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))

#
TE_CDS %>% group_by(Category) %>% summarise(median=median(RNA)/7.82)
TE_CDS %>% group_by(Category) %>% summarise(median=median(TE)/1.19)
half_lives %>% group_by(Category) %>% summarise(median=median(Mean)/1.90)
shoot_proteomics2 %>% group_by(Category) %>% summarise(median=median(norm_MS_count)/0.0468)
#######################
pRNA <- ggplot(TE_CDS, aes(x=RNA,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,40)+xlab("RNA (TPM)")+ylab("Fraction of Genes")+labs(tag = "A") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12)) + 
  ggtitle("RNA")+
  scale_colour_discrete("Category",labels=c("TuORFs","UuORFs","Others"))
pRNA 

pTE <- ggplot(TE_CDS, aes(x=TE,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,4)+xlab("TE")+ylab("Fraction of Genes")+labs(tag = "B") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12)) + 
  ggtitle("TE")+
  scale_colour_discrete("Category",labels=c("TuORFs","UuORFs","Others"))
pTE

pHL <- ggplot(half_lives, aes(x=Mean,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,7)+xlab("Half-lives (Hours)")+ylab("Fraction of Genes")+labs(tag = "C") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12)) + 
  ggtitle("Half-lives") +
  scale_colour_discrete("Category",labels=c("TuORFs","UuORFs","Others"))
pHL

pProtein <- ggplot(shoot_proteomics2, aes(x=norm_MS_count,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.5)+xlab("Protein abundance")+ylab("Fraction of Genes")+labs(tag = "D") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12),
                          legend.title = element_text(size=12)) + 
  ggtitle("Protein abundance") +
  scale_colour_discrete("Category",labels=c("TuORFs","UuORFs","Others"))
pProtein

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pProtein)

p3 =grid.arrange(pRNA + theme(legend.position="none"),
                 pTE + theme(legend.position="none"),
                 pHL + theme(legend.position="none"),
                 pProtein + theme(legend.position="none"),
                 p2_legend,
                 layout_matrix=rbind(c(1,2,3,4,5)),
                 top=textGrob("",gp = gpar(fontsize = 14, fontface = "bold")))

ggsave("~/Desktop/atRTD3/Figure2A-D ECDF.pdf",
       plot = p3,
       units = "in",
       width = 13,
       height = 3.5)
