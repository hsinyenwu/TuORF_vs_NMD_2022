#############################################################################################
# Figure S8: Length distribution of CDS for Ribo-seq uORFs, predicted uORFs and no uORF genes
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

fiveUTRByTx_length <- sum(width(fiveUTRByTx))
fUTRSbyGene_length_df <- data.frame(gene_id=substr(names(fiveUTRByTx_length),1,9),length_nt=unname(fiveUTRByTx_length))
fUTRSbyGene_length_df <- fUTRSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

threeUTRByTx_length <- sum(width(threeUTRByTx))
tUTRSbyGene_length_df <- data.frame(gene_id=substr(names(threeUTRByTx_length),1,9),length_nt=unname(threeUTRByTx_length))
tUTRSbyGene_length_df <- fUTRSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

CDSbyGene <- cdsBy(txdb,by="gene")
CDSbyGene_length <- sum(width(CDSbyGene))
CDSbyGene_length_df <- data.frame(gene_id=names(CDSbyGene_length),length_nt=unname(CDSbyGene_length))
CDSbyGene_length_df <- CDSbyGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

CDSbyTx <- cdsBy(txdb,by="tx",use.names=T)

# exonByGene <- exonsBy(txdb,by='gene')
exonByGene_length <- sum(width(exonByGene))
exonByGene_length_df <- data.frame(gene_id=names(exonByGene_length),length_nt=unname(exonByGene_length))
exonbyGene_length_df <- exonByGene_length_df %>% filter(gene_id %in% fiveUTRByGeneNames)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")

Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id)]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id
#Gene without uORfs but with 5'UTR
Genes_wo_uORFs <- fiveUTRByGeneNames[!(fiveUTRByGeneNames %in% c(ORF_max_filt_uORF$gene_id,Genes_w_only_seq_predicted_uORFs))]
head(Genes_wo_uORFs)
ORF_max_filt_predicted_uORF <- ORF_max_filt %>% filter(gene_id %in% Genes_w_only_seq_predicted_uORFs)
ORF_max_filt_no_uORF <- ORF_max_filt %>% filter(gene_id %in% Genes_wo_uORFs)

fUTRSbyGene_length_df$Category <- ifelse(fUTRSbyGene_length_df$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"TuORFs",ifelse(fUTRSbyGene_length_df$gene_id %in% Genes_w_only_seq_predicted_uORFs, "UuORFs", "Others"))
fUTRSbyGene_length_df$Category <- factor(fUTRSbyGene_length_df$Category, levels=c("TuORFs", "UuORFs", "Others"), labels=c("TuORFs", "UuORFs", "Others"))
table(fUTRSbyGene_length_df$Category)

tUTRSbyGene_length_df$Category <- ifelse(tUTRSbyGene_length_df$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"TuORFs",ifelse(tUTRSbyGene_length_df$gene_id %in% Genes_w_only_seq_predicted_uORFs, "UuORFs", "Others"))
tUTRSbyGene_length_df$Category <- factor(tUTRSbyGene_length_df$Category, levels=c("TuORFs", "UuORFs", "Others"), labels=c("TuORFs", "UuORFs", "Others"))
table(tUTRSbyGene_length_df$Category)

CDSbyGene_length_df$Category <- ifelse(CDSbyGene_length_df$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"TuORFs",ifelse(CDSbyGene_length_df$gene_id %in% Genes_w_only_seq_predicted_uORFs, "UuORFs", "Others"))
CDSbyGene_length_df$Category <- factor(CDSbyGene_length_df$Category, levels=c("TuORFs", "UuORFs", "Others"), labels=c("TuORFs", "UuORFs", "Others"))
table(CDSbyGene_length_df$Category)

exonByGene_length_df$Category <- ifelse(exonByGene_length_df$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"TuORFs",ifelse(exonByGene_length_df$gene_id %in% Genes_w_only_seq_predicted_uORFs, "UuORFs", "Others"))
exonByGene_length_df$Category <- factor(exonByGene_length_df$Category, levels=c("TuORFs", "UuORFs", "Others"), labels=c("TuORFs", "UuORFs", "Others"))
table(exonByGene_length_df$Category)


#
fUTRSbyGene_length_df %>% group_by(Category) %>% summarise(median=median(length_nt)/126)
CDSbyGene_length_df %>% group_by(Category) %>% summarise(median=median(length_nt)/1071)
tUTRSbyGene_length_df %>% group_by(Category) %>% summarise(median=median(length_nt)/244)
exonByGene_length_df %>% group_by(Category) %>% summarise(median=median(length_nt)/1367)
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
                 ylim=c(-1,900),
                 legend = "right")+
  rotate_x_text(angle = 45) +
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

vsteps <- yvalue(ystart=820,stepsize=40,num=3)

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
                 ylim=c(-1,4100),
                 legend = "right")+
  rotate_x_text(angle = 45) +
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

vsteps <- yvalue(ystart=3800,stepsize=150,num=3)

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
                 ylim=c(-1,700),
                 legend = "right")+
  rotate_x_text(angle = 45) +
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

vsteps <- yvalue(ystart=620,stepsize=40,num=3)

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
                 ylim=c(-1,5000),
                 legend = "right")+
  rotate_x_text(angle = 45) +
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

vsteps <- yvalue(ystart=4600,stepsize=200,num=3)

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
                  top=textGrob("FigureS8: Length 5UTR 3UTR CDS uORF genes", gp=gpar(fontsize=12), x = 0, hjust = 0))

ggsave("~/Desktop/Figure S8 Length 5UTR 3UTR CDS uORF genes.pdf",
       plot = p3,
       units = "in",
       width = 10.5,
       height = 4)

