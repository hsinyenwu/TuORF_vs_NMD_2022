#############################################################################################
# Figure S3 Root mass-spec protein quantification for TuORF vs UuORF vs Other genes
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

# RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
# Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
# TE <- Ribo$TPM/RNA$TPM
#
# TE_CDS <- data.frame(gene_id=RNA$gene_id,
#                      RNA=RNA$TPM,
#                      Ribo=Ribo$TPM,
#                      TE=TE)
# #genes with 5'UTR
fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,9)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)
# write.table(fiveUTRByGeneNames_df,"~/Desktop/uORFs_miRNA/genes_w_5'UTR.tsv",quote = F, row.names = FALSE,col.names = TRUE, sep="\t")
# TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)
# # Remove not translated transcripts
# TE_CDS <- TE_CDS %>% filter(TE>0)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id
#
# TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(TE_CDS$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others"))
# TE_CDS$Category <- factor(TE_CDS$Category, levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))


# Quantitative peoteomics for root
root_proteomics <- read.xlsx("~/Desktop/uORFs_miRNA/Justin_Proteomics2017/data_excel/pmic12946-sup-0005-tables9.xlsx",sheet = 3)
root_proteomics$gene_id <- substr(root_proteomics$Protein.IDs,1,9)
root_proteomics2 <- root_proteomics %>% dplyr::select(gene_id,`MS/MS.count`,`Sequence.length`)
colnames(root_proteomics2) <-c("gene_id","MS_count","Pept_length")
root_proteomics2 <- root_proteomics2 %>% filter(gene_id %in% fiveUTRByGeneNames)
root_proteomics2$norm_MS_count <- root_proteomics2$MS_count/root_proteomics2$Pept_length
# root_proteomics2$Category <- ifelse(root_proteomics2$gene_id %in% CPuORFs$gene_id, "CPuORFs",ifelse(root_proteomics2$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(root_proteomics2$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others")))
root_proteomics2$Category <- ifelse(root_proteomics2$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq",ifelse(root_proteomics2$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted", "Others"))
root_proteomics2$Category <- factor(root_proteomics2$Category,levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))

#
root_proteomics2 %>% group_by(Category) %>% summarise(median=median(norm_MS_count)/0.0358)

#######################
# protein abundance boxplot
#######################
root_proteomics2 <- root_proteomics2[!is.infinite(root_proteomics2$norm_MS_count),]

# default p-value adjustement method: Holm
stat.test <- root_proteomics2 %>%
  wilcox_test(norm_MS_count ~Category) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(root_proteomics2,
                 x = "Category", y = "norm_MS_count",
                 color = "Category",
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.01,0.36),
                 legend = "right",
                 ylab="Normalized MS count") +
  geom_hline(yintercept=median(root_proteomics2$norm_MS_count[root_proteomics2$Category=="Predicted"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) +
  labs(tag = "D") +
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_color_manual(values = c("Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF"))
#legend key font size

bxp

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

yvalue <- function(ystart,stepsize,num) {
  Ev <- c()
  for(i in 1:num){
    Ev[i] <- ystart + stepsize*(i-1)
  }
  return(Ev)
}

vsteps <- yvalue(ystart=0.26,stepsize=0.02,num=3)

pProtein_box <- bxp +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.01,size = 2.5)
pProtein_box

