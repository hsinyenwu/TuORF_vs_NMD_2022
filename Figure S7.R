```
#############################################################################################
# Figure 3: Length distribution of CDS for Ribo-seq uORFs, predicted uORFs and no uORF genes
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
threeUTRByTx <- threeUTRsByTranscript(txdb,use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 22136
fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,9)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=T,groupByTx=F, minimumLength=9)
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])))
fiveUTR_ORFs_gene_names <- substr(fiveUTR_ORFs_tx_names,1,9)

threeUTRByTx_length <- sum(width(threeUTRByTx))
tUTRSbyGene_length_df <- data.frame(gene_id=substr(names(threeUTRByTx_length),1,9),length_nt=unname(threeUTRByTx_length))
tUTRSbyGene_length_df <- fUTRSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

fiveUTRByTx_length <- sum(width(fiveUTRByTx))
fUTRSbyGene_length_df <- data.frame(gene_id=substr(names(fiveUTRByTx_length),1,9),length_nt=unname(fiveUTRByTx_length))
fUTRSbyGene_length_df <- fUTRSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

CDSbyGene <- cdsBy(txdb,by="gene")
CDSbyGene_length <- sum(width(CDSbyGene))
CDSbyGene_length_df <- data.frame(gene_id=names(CDSbyGene_length),length_nt=unname(CDSbyGene_length))
CDSbyGene_length_df <- CDSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

CDSbyTx <- cdsBy(txdb,by="tx",use.names=T)

# exonByGene <- exonsBy(txdb,by='gene')
exonByGene_length <- sum(width(exonByGene))
exonByGene_length_df <- data.frame(gene_id=names(exonByGene_length),length_nt=unname(exonByGene_length))
exonbyGene_length_df <- exonByGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

# #RiboTaper output
# ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
# table(ORF_max_filt$category)
# ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
# 
# NMD_target_w_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id)]
# NMD_target_wo_uORFs <- ORF_max_filt_uORF$gene_id
# #Gene without uORfs but with 5'UTR
# Genes_wo_uORFs <- fiveUTRByGeneNames[!(fiveUTRByGeneNames %in% c(ORF_max_filt_uORF$gene_id,NMD_target_w_uORFs))]
# head(Genes_wo_uORFs)
# ORF_max_filt_predicted_uORF <- ORF_max_filt %>% filter(gene_id %in% NMD_target_w_uORFs)
# ORF_max_filt_no_uORF <- ORF_max_filt %>% filter(gene_id %in% Genes_wo_uORFs)
########################################################################
# Figure 4. NMD targets, uORF genes, other genes RNA TE Half-lives Protein
########################################################################
#

#Load 5UTR containing gene list

fiveUTRByGeneNames <- read.delim("~/Desktop/uORFs_miRNA/genes_w_5'UTR.tsv",header=T)
fiveUTRByGeneNames <- fiveUTRByGeneNames$gene_id

#Load RNA-seq Ribo-seq datasets from RSEM output
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

#Load RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
ORF_max_filt_1uORF <- ORF_max_filt_uORF %>% group_by(gene_id) %>% filter(n()==1)
nrow(ORF_max_filt_1uORF) #[1] 1515
ORF_max_filt_1uORF$pept_length <- nchar(ORF_max_filt_1uORF$ORF_pept)
sum(ORF_max_filt_1uORF$pept_length>=35) #308
sum(ORF_max_filt_1uORF$pept_length<35) #1207
ORF_max_filt_1uORF_gt_35 <- ORF_max_filt_1uORF %>% filter(pept_length>=35)
dim(ORF_max_filt_1uORF_gt_35)
ORF_max_filt_1uORF_gt_35_NMD <- ORF_max_filt_1uORF_gt_35 %>% filter(gene_id %in% NMD$gene_id)
dim(ORF_max_filt_1uORF_gt_35_NMD)

#Load NMD targets and classify
NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
# NMD <- read.xlsx("~/Desktop/uORFs_miRNA/Lee_eta_al_EJC_NMD_candidates.xlsx")
colnames(NMD) <- "gene_id"

NMD_target_wo_uORFs <- setdiff(NMD$gene_id,ORF_max_filt_uORF$gene_id)
length(NMD_target_wo_uORFs) #[1] 284
NMD_target_w_uORFs <- intersect(NMD$gene_id,ORF_max_filt_uORF$gene_id)
length(NMD_target_w_uORFs) #[1] 49
uORFs_not_NMD_targets <- setdiff(ORF_max_filt_uORF$gene_id,NMD$gene_id)
length(uORFs_not_NMD_targets) #[1] 1719
Others <- setdiff(ORF_max_filt$gene_id,c(NMD$gene_id,ORF_max_filt_uORF$gene_id))
length(Others) #19461

ORF_max_filt_1uORF <- ORF_max_filt_uORF %>% group_by(gene_id) %>% filter(n()==1)
ORF_max_filt_1uORF_NMD <- ORF_max_filt_1uORF %>% filter(gene_id %in% NMD_target_w_uORFs)
nrow(ORF_max_filt_1uORF_NMD) #37
ORF_max_filt_1uORF_NMD$pept_length <- nchar(ORF_max_filt_1uORF_NMD$ORF_pept)
sum(ORF_max_filt_1uORF_NMD$pept_length>=35) #15
sum(ORF_max_filt_1uORF_NMD$pept_length<35) #22

TE_CDS_NMD <- TE_CDS %>% filter(gene_id %in% NMD$gene_id)
dim(TE_CDS_NMD)
# library(clipr)
# write_clip(TE_CDS_NMD$gene_id)

# Expression levels (for plotting RNA-seq and TE)
TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(TE_CDS$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(TE_CDS$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(TE_CDS$Category)

TE_CDS_atuORF_NMDt <- TE_CDS %>% filter(Category=="NMD_target_w_uORFs")
library(clipr)
write_clip(TE_CDS_atuORF_NMDt$gene_id)


fUTRSbyGene_length_df$Category <- ifelse(fUTRSbyGene_length_df$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(fUTRSbyGene_length_df$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(fUTRSbyGene_length_df$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
fUTRSbyGene_length_df$Category <- factor(fUTRSbyGene_length_df$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(fUTRSbyGene_length_df$Category)

tUTRSbyGene_length_df$Category <- ifelse(tUTRSbyGene_length_df$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(tUTRSbyGene_length_df$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(tUTRSbyGene_length_df$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
tUTRSbyGene_length_df$Category <- factor(tUTRSbyGene_length_df$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(tUTRSbyGene_length_df$Category)

CDSbyGene_length_df$Category <- ifelse(CDSbyGene_length_df$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(CDSbyGene_length_df$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(CDSbyGene_length_df$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
CDSbyGene_length_df$Category <- factor(CDSbyGene_length_df$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(CDSbyGene_length_df$Category)

exonByGene_length_df$Category <- ifelse(exonByGene_length_df$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(exonByGene_length_df$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(exonByGene_length_df$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
exonByGene_length_df$Category <- factor(exonByGene_length_df$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(exonByGene_length_df$Category)

#
fUTRSbyGene_length_df %>% group_by(Category) %>% dplyr::summarise(Median=median(length_nt))

fUTRSbyGene_length_df %>% group_by(Category) %>% dplyr::summarise(median=median(length_nt)/162)
CDSbyGene_length_df %>% group_by(Category) %>% dplyr::summarise(median=median(length_nt)/1080)
tUTRSbyGene_length_df %>% group_by(Category) %>% dplyr::summarise(median=median(length_nt)/248)
exonByGene_length_df %>% group_by(Category) %>% dplyr::summarise(median=median(length_nt)/1458)
#######################
####   Box plot    ####
#######################

#######################
# 5'UTR length boxplot
#######################

stat.test <- fUTRSbyGene_length_df %>%
  wilcox_test(length_nt ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(fUTRSbyGene_length_df, 
                 x = "Category", y = "length_nt", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,1300),
                 legend = "right")+
  rotate_x_text(angle = 45) + 
  geom_hline(yintercept=median(fUTRSbyGene_length_df$length_nt[fUTRSbyGene_length_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  labs(tag = "A") + 
  ylab("5' UTR length (nt)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) 
# bxp

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

vsteps <- yvalue(ystart=1020,stepsize=50,num=6)

p5UTR_length_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
p5UTR_length_box

#######################
# CDS length boxplot
#######################

stat.test <- CDSbyGene_length_df %>%
  wilcox_test(length_nt ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(CDSbyGene_length_df, 
                 x = "Category", y = "length_nt", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,4500),
                 legend = "right")+
  rotate_x_text(angle = 45) + 
  geom_hline(yintercept=median(CDSbyGene_length_df$length_nt[CDSbyGene_length_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  labs(tag = "B") + 
  ylab("CDS length (nt)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) 
# bxp

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

vsteps <- yvalue(ystart=3650,stepsize=160,num=6)

pCDS_length_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pCDS_length_box

#######################
# 3'UTR length boxplot
#######################

stat.test <- tUTRSbyGene_length_df %>%
  wilcox_test(length_nt ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(tUTRSbyGene_length_df, 
                 x = "Category", y = "length_nt", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,750),
                 legend = "right")+
  rotate_x_text(angle = 45) + 
  geom_hline(yintercept=median(tUTRSbyGene_length_df$length_nt[tUTRSbyGene_length_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  labs(tag = "C") + 
  ylab("3' UTR length (nt)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) 
# bxp

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

vsteps <- yvalue(ystart=580,stepsize=30,num=6)

p3UTR_length_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
p3UTR_length_box

#######################
# transcript (exon) length boxplot
#######################

stat.test <- exonByGene_length_df %>%
  wilcox_test(length_nt ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(exonByGene_length_df, 
                 x = "Category", y = "length_nt", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,5800),
                 legend = "right")+
  rotate_x_text(angle = 45) + 
  geom_hline(yintercept=median(exonByGene_length_df$length_nt[exonByGene_length_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  labs(tag = "D") + 
  ylab("Transcript length (nt)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) 
# bxp

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

vsteps <- yvalue(ystart=4430,stepsize=250,num=6)

pTx_length_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pTx_length_box

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pTx_length_box)

p3 = grid.arrange(p5UTR_length_box + theme(legend.position="none"),
                  pCDS_length_box + theme(legend.position="none"),
                  p3UTR_length_box + theme(legend.position="none"),
                  pTx_length_box+ theme(legend.position="none"),
                  p2_legend,
                  layout_matrix=rbind(c(1,2,3,4,5)),
                  top=textGrob("Figure X: Length 5UTR 3UTR CDS NMD targets and uORF genes", gp=gpar(fontsize=12), x = 0, hjust = 0))

ggsave("~/Desktop/NMD paper/Figure Length 5UTR 3UTR CDS NMD targets and uORF genes.pdf",
       plot = p3,
       units = "in",
       width = 12,
       height = 5)
```
