##################################################
# Figure S6: CPuORFs vs Other genes              #                         
##################################################
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
# txdb <- makeTxDbFromGFF("~/Desktop/uORFs_miRNA/Araport11_20181206_max_isoform.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206_expressed.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")

exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 37994

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=F,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]][1])))
length(fiveUTR_ORFs_tx_names) #[1] 17125
fiveUTR_ORFs_gene_names <- substr(fiveUTR_ORFs_tx_names,1,9)

RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
TE <- Ribo$TPM/RNA$TPM

# CPuORFs <- read.xlsx("~/Desktop/uORFs_miRNA/TableS2_uORFs_on_CPuORF_containg_5UTR_comprehensive_sing_CPuORF_genes_only_colored.xlsx")
CPuORFs <- read.xlsx("~/Desktop/CTRL_v1/CPuORF_collection_3_manully_fix_b.xlsx",sheet =5)

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
TE_CDS <- TE_CDS %>% filter(RNA>0,TE>0)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093 
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% CPuORFs$gene_id,"CPuORFs",ifelse(TE_CDS$gene_id %in% fiveUTR_ORFs_gene_names, "Other uORFs", "Others"))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("CPuORFs", "Other uORFs", "Others"), labels=c("CPuORFs", "Other uORFs","Others"))

TE_CDS %>% group_by(Category) %>% summarise(median=median(RNA)/7.66)

# TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% CPuORFs$gene_id,"CPuORFs",ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs, "TuORFs", "Others"))
# table(TE_CDS$Category)
# # TuORFs CPuORFs  Others 
# # 1693      84   19859 
# TE_CDS$Category <- factor(TE_CDS$Category, levels=c("CPuORFs", "TuORFs", "Others"), labels=c("CPuORFs", "TuORFs","Others"))

#######################
####   Box plot    ####
#######################

#######################
# RNA boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(log2(RNA) ~ Category, data = TE_CDS)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(log2(RNA) ~ Category, data = TE_CDS, var.equal = FALSE)
# p-value < 2.2e-16

# default p-value adjustement method: Holm
stat.test <- TE_CDS %>%
  wilcox_test(RNA ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

TE_CDS %>% group_by(Category) %>% summarise(Median=median(RNA))
# Create a box plot
bxp <- ggboxplot(TE_CDS, 
                 x = "Category", y = "RNA", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,82),
                 legend = "right") + 
  geom_hline(yintercept=median(TE_CDS$RNA[TE_CDS$Category=="Predicted"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "A") + 
  ylab("RNA (TPM)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_color_manual(values = c("CPuORFs" = "#F8766D", "Other uORFs"="#00BA38", "Others"="#619CFF")) 

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

vsteps <- yvalue(ystart=72,stepsize=4,num=3)


pRNA_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0001, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pRNA_box
