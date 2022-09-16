####################################################################################################################
# Figure 2 A-D Cumulative plot (ECDF) for TuORF, UuORF and other genes (i.e., no potential uORFs but with annotated 5'UTR 
# 2A steady-state mRNA levels
# 2B TE (translation efficiency)
# 2C mRNA half-lives
# 2D protein abundance
# 2E-H Same as 2A-D but present with boxplots
# 2I 
####################################################################################################################

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
TE_CDS %>% group_by(Category) %>% dplyr::summarise(median=median(RNA)/7.82)
TE_CDS %>% group_by(Category) %>% dplyr::summarise(median=median(TE)/1.19)
half_lives %>% group_by(Category) %>% dplyr::summarise(median=median(Mean)/1.90)
shoot_proteomics2 %>% group_by(Category) %>% dplyr::summarise(median=median(norm_MS_count)/0.0468)

#########################
# Plot Figure 2A-D ECDF #
#########################
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

#ggsave("~/Desktop/atRTD3/Figure2A-D ECDF.pdf",
#       plot = p3,
#       units = "in",
#       width = 13,
#       height = 3.5)

#######################
#      Box plot       #
#     Figure 2E-H     #
#######################

#######################
# Figure 2E                      
# RNA boxplot
#######################

##############################
# default p-value adjustement method: Holm
stat.test <- TE_CDS %>%
  wilcox_test(RNA ~Category) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(TE_CDS,
                 x = "Category", y = "RNA",
                 color = "Category",
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,60),
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
  scale_color_manual(values = c("Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF"))

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

vsteps <- yvalue(ystart=51,stepsize=3,num=3)


pRNA_box <- bxp +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0001, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pRNA_box

#######################
# Figure 2F
# TE boxplot
#######################

TE_CDS %>% group_by(Category) %>% dplyr::summarise(Median=median(TE))
TE_CDS <- TE_CDS[!is.infinite(TE_CDS$TE),]

stat.test <- TE_CDS %>%
  wilcox_test(TE ~ Category) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(TE_CDS,
                 x = "Category", y = "TE",
                 color = "Category",
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.1,3.7),
                 legend = "right") +
  geom_hline(yintercept=median(TE_CDS$TE[TE_CDS$Category=="Predicted"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) +
  labs(tag = "B") +
  ylab("TE (TPM)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
  scale_color_manual(values = c("Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF"))

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

vsteps <- yvalue(ystart=2.5,stepsize=0.25,num=3)

pTE_box <- bxp +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0001, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.1,size = 2.5)
pTE_box

#######################
# Figure 2G
# Half-life boxplot
#######################
half_lives <- half_lives[!is.infinite(half_lives$Mean),]

# default p-value adjustement method: Holm
stat.test <- half_lives %>%
  wilcox_test(Mean ~Category) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(half_lives,
                 x = "Category", y = "Mean",
                 color = "Category",
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.1,6),
                 legend = "right",
                 ylab="Mean half-life") +
  geom_hline(yintercept=median(half_lives$Mean[half_lives$Category=="Predicted"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) +
  labs(tag = "C") +
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))+
  scale_color_manual(values = c("Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF"))  #legend key font size

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

vsteps <- yvalue(ystart=4.3,stepsize=0.35,num=3)

pHL_box <- bxp +
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.005, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.1,size = 2.5)
pHL_box

#############################
# Figure 2H                 #
# protein abundance boxplot #
#############################
shoot_proteomics2 <- shoot_proteomics2[!is.infinite(shoot_proteomics2$norm_MS_count),]

# default p-value adjustement method: Holm
stat.test <- shoot_proteomics2 %>%
  wilcox_test(norm_MS_count ~Category) %>%
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(shoot_proteomics2,
                 x = "Category", y = "norm_MS_count",
                 color = "Category",
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.01,0.36),
                 legend = "right",
                 ylab="Normalized MS count") +
  geom_hline(yintercept=median(shoot_proteomics2$norm_MS_count[shoot_proteomics2$Category=="Predicted"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
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

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pProtein_box)

p3 <- grid.arrange(pRNA_box + theme(legend.position="none"),
                   pTE_box + theme(legend.position="none"),
                   pHL_box + theme(legend.position="none"),
                   pProtein_box+ theme(legend.position="none"),
                   p2_legend,
                   layout_matrix=rbind(c(1,2,3,4,5)),
                   top=textGrob("Figure 2.Riboseq identified uORF genes vs predicted uORF genes vs other genes", gp=gpar(fontsize=13), x = 0, hjust = 0))
p3

# ggsave("~/Desktop/atRTD3/Figure2E-H_boxplot.pd",
#        plot = p3,
#        units = "in",
#        width = 13,
#        height = 4)


###################################
# 2I Boxplot for multiple tissue/stages
###################################
# dataset 6: Arabidopsis tissue atlas E-MTAB-7978
dataset6 <- read.xlsx("~/Desktop/uORFs_miRNA/Figures_May052021/E-MTAB-7978-query-results.tpms.xlsx",sheet = 1,startRow =5)
head(dataset6)
library(tidyr)
dataset6 <- dataset6 %>% replace(is.na(.), 0)
colnames(dataset6)[1] <- "gene_id"
dim(dataset6) #[1] 20897     7
dataset6$TPM_MEAN <- rowMeans(dataset6[,3:ncol(dataset6)])
head(dataset6,20)
dataset6 <- dataset6 %>% filter(gene_id %in% fiveUTRByGeneNames)
dataset6$Category <- ifelse(dataset6$gene_id %in% ORF_max_filt_uORF_id,"TuORFs",ifelse(dataset6$gene_id %in% seq_defined_uORFs, "UuORFs", "Others"))
dataset6$Category <- factor(dataset6$Category,levels=c("TuORFs", "UuORFs", "Others"), labels=c("TuORFs", "UuORFs","Other genes"))

table(dataset6$Category)

dataset6_mean_by_group <- aggregate(dataset6[,3:(ncol(dataset6)-1)], list(dataset6$Category), FUN=median) %>% as.data.frame()
head(dataset6_mean_by_group)
colnames(dataset6_mean_by_group)[1] <- "Category"

dataset6_mean_by_group$TPM_MEAN <- rowMeans(dataset6_mean_by_group[,3:ncol(dataset6_mean_by_group)])

library(reshape2)
dataset6_mean_by_group_melt <- melt(dataset6_mean_by_group)

dataset6m <- melt(dataset6)
p<- ggplot(dataset6m, aes(x=variable, y=value, color=Category)) + ylim(0,85) +
  stat_boxplot(geom ='errorbar') +
  geom_boxplot(outlier.shape = NA) + theme_classic() +
  xlab("Samples") + ylab("RNA (TPM)")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p

#ggsave("~/Desktop/Figure 2I.pdf",
       # plot = p,
       # units = "in",
       # width = 13,
       # height = 6)                      
                      
                     
