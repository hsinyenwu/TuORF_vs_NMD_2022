####################################################################################################################
# Figure 2 A-D Cumulative plot for TuORF, UuORF and other genes (i.e., no potential uORFs but with annotated 5'UTR 
# 2A steady-state mRNA levels
# 2B TE (translation efficiency)
# 2C mRNA half-lives
# 2D protein abundance
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

ggsave("~/Desktop/Figure2A-D ECDF.pdf",
       plot = p3,
       units = "in",
       width = 13,
       height = 3.5)



                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
                      
##################################################
# Figure 2: RT-uORFs vs Sequence predicted uORFs #
# RNA/Ribo/TE/HL/Protein levels                  #
# Boxplot and U test                             #                              
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

#######################
####   Box plot    ####
#######################

#######################
# RNA boxplot
#######################

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
# TE boxplot
#######################

TE_CDS %>% group_by(Category) %>% summarise(Median=median(TE))
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
  scale_color_manual(values = c("CPuORFs"="purple","Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF"))  #legend key font size

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

#######################
# protein abundance boxplot
#######################
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
  scale_color_manual(values = c("CPuORFs"="purple","Riboseq" = "#F8766D", "Predicted"="#00BA38", "Others"="#619CFF")) 
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

# ggsave("~/Desktop/uORFs_miRNA/Figures_May052021/final/Figure2_Riboseq_vs_seq_uORFs_5b.pdf",
#        plot = p3,
#        units = "in",
#        width = 13,
#        height = 4)


################
# Figure 2E-F
################
#rm(list=ls())
library(dplyr)
library(GenomicRanges)
library(GenomicFeatures)
library(Biostrings)
library(Rsamtools)
library(seqLogo)

ORFs_max_filt <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/ORFs_max_filt_20181118",header=T,sep="\t",stringsAsFactors = F,quote = "")
nrow(ORFs_max_filt)
ORFs_max_filt_uORFs <- ORFs_max_filt %>% filter(category=="uORF")
nrow(ORFs_max_filt_uORFs)
ORFs_max_filt_w_5UTR <- ORFs_max_filt %>% filter(annotated_start >1 & category =="ORFs_ccds")
dim(ORFs_max_filt_w_5UTR) #[1] 14732  38
table(ORFs_max_filt$category)

# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 163        354      20659          7       1329 


#Find genes with high TE and check the translation start sites
# RNAseq and riboseq TPM and TE
R1RNA <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/RNA.genes.results/R1.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
R2RNA <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/RNA.genes.results/R2.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
R3RNA <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/RNA.genes.results/R3.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
head(R1RNA)
dim(R1RNA) #[1] 36086 7
RNA_TPM <- data.frame(gene_id=R1RNA$gene_id,R1=R1RNA$TPM,R2=R2RNA$TPM,R3=R3RNA$TPM)
RNA_TPM$TPM_ave <- round((RNA_TPM$R1+RNA_TPM$R2+RNA_TPM$R3)/3,2)
head(RNA_TPM)

R1ribo <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/ribo.genes.results/R1.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
R2ribo <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/ribo.genes.results/R2.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
R3ribo <- read.delim("~/Desktop/Tomato_paper_version/Tomato_Root_20181118/ribo.genes.results/R3.genes.results",header=T,sep="\t",stringsAsFactors = F,quote = "")
head(R1ribo)
dim(R1ribo) #[1] 36089 7
ribo_TPM <- data.frame(gene_id=R1ribo$gene_id,R1=R1ribo$TPM,R2=R2ribo$TPM,R3=R3ribo$TPM)
ribo_TPM$TPM_ave <- round((ribo_TPM$R1+ribo_TPM$R2+ribo_TPM$R3)/3,2)
head(ribo_TPM)
# plot(log2(RNA_TPM$TPM_ave),log2(ribo_TPM$TPM_ave),xlab="log2(RNA_TPM)",ylab="log2(Ribo_TPM)")
TE <- ribo_TPM$TPM_ave/RNA_TPM$TPM_ave
TE_CDS <- data.frame(gene_id=R1ribo$gene_id,RNA_TPM_ave=RNA_TPM$TPM_ave,ribo_TPM_ave=ribo_TPM$TPM_ave,TE=TE)

TE_CDS <- TE_CDS %>% filter(RNA_TPM_ave>0 & TE>0)
colnames(TE_CDS) <- c("gene_id","RNA","Ribo","TE")

# TE_CDS <- TE_CDS[TE_CDS$gene_id %in% ORFs_max_filt_w_5UTR$gene_id,]

# colnames(ORFs_max_filt)
FA <- FaFile("~/Desktop/Tomato/Reference/seqs_srcdir/S_lycopersicum_3.00.fa")
txdb <- makeTxDbFromGFF("~/Desktop/Tomato_paper_version//Tomato_Root_20181118/SupTables/Tomato_Root_ITAG_vs_de_novo_10312018.ixyous+ITAG3.2_rm_start_stop.gtf",format="gtf", dataSource="ITAG3.2",organism="Solanum lycopersicum")
exonByGene <- exonsBy(txdb,by='gene')
exonByTx <- exonsBy(txdb,by='tx',use.names=T)
fiveUTRByTx <- fiveUTRsByTranscript(txdb,use.names=T)
fiveUTR_seqs <- extractTranscriptSeqs(FA,fiveUTRByTx)
length(fiveUTR_seqs) #[1] 17868

#genes with 5'UTR
fiveUTRByGeneNames <- substr(names(fiveUTRByTx),1,16)
fiveUTRByGeneNames_df <- data.frame(gene_id=fiveUTRByGeneNames)

fiveUTR_ORFs <- findMapORFs(fiveUTRByTx, fiveUTR_seqs,startCodon = "ATG",longestORF=F,groupByTx=F, minimumLength=9)
length(unique(fiveUTR_ORFs)) #42328
fiveUTR_ORFs_tx_names <- unlist(lapply(1:length(fiveUTR_ORFs),function(x) names(fiveUTR_ORFs[[x]])[1]))
length(fiveUTR_ORFs_tx_names) #42328
length(unique(fiveUTR_ORFs_tx_names))#[1] 9166
fiveUTR_ORFs_gene_names <- substr(unique(fiveUTR_ORFs_tx_names),1,16)
length(fiveUTR_ORFs_gene_names) #9166

tail(fiveUTR_ORFs_gene_names)
sum(grepl("Solyc",fiveUTR_ORFs_gene_names)) #9166

Only_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names %in% ORFs_max_filt_uORFs$gene_id)]

TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% ORFs_max_filt_uORFs$gene_id,"Riboseq",ifelse(TE_CDS$gene_id %in% Only_predicted_uORFs, "Predicted", "Others"))
table(TE_CDS$Category)
# Others Predicted   Riboseq 
# 7069      5995      1124 
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("Riboseq", "Predicted", "Others"), labels=c("Riboseq", "Predicted","Others"))

#
TE_CDS %>% group_by(Category) %>% summarise(median=median(RNA)/16.5)
TE_CDS %>% group_by(Category) %>% summarise(median=median(TE)/0.955)

#######################
####   Box plot    ####
#######################
#######################
#######################
#######################
library(ggpubr)
library(rstatix)
library(gridExtra)
library(grid)
library(EnvStats)
#######################
# RNA boxplot
#######################

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
                 ylim=c(-1,100),
                 legend = "right") + 
  geom_hline(yintercept=median(TE_CDS$RNA[TE_CDS$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "A") + 
  ylab("RNA (TPM)")+
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

vsteps <- yvalue(ystart=88,stepsize=4,num=3)

# # Add p-values onto the box plots
# stat.test <- stat.test %>%
#   add_xy_position(x = "Category", dodge = 0.8)

pRNA_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.00002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pRNA_box

#######################
# TE boxplot
#######################

TE_CDS <- TE_CDS[!is.infinite(TE_CDS$TE),]
# TE_CDS <-  TE_CDS[TE_CDS$TE<5,]
# default p-value adjustement method: Holm
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
                 ylim=c(-0.1,3.8),
                 legend = "right") + 
  geom_hline(yintercept=median(TE_CDS$TE[TE_CDS$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "B") + 
  ylab("TE (TPM)")+
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

vsteps <- yvalue(ystart=2.3,stepsize=0.25,num=3)
# vsteps <- yvalue(ystart=3,stepsize=0.25,num=6)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

pTE_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.00002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.1,size = 2.5)
pTE_box


get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pTE_box)

p3= grid.arrange(pRNA_box + theme(legend.position="none"),
                 pTE_box + theme(legend.position="none"),
                 p2_legend,
                 layout_matrix=rbind(c(1,2,3)),
                 top=textGrob("Figure 2. Riboseq validated vs sequence predicted uORFs vs others in tomato", gp=gpar(fontsize=13)))


# ggsave("~/Desktop/uORFs_miRNA/Figures_May052021/final/Figure2_tomato_RNA_TE_box.pdf",
#        plot = p3,
#        units = "in",
#        width = 7,
#        height = 3)


#############################################################################################
# Figure 2H: Expression levels for ATuORFs, PTuORFs and no uORF genes
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

# LOad RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")

Genes_w_only_seq_predicted_uORFs <- fiveUTR_ORFs_gene_names[!(fiveUTR_ORFs_gene_names%in%ORF_max_filt_uORF$gene_id)]
Genes_w_only_RiboTaper_defined_uORFs <- ORF_max_filt_uORF$gene_id

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

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% Genes_w_only_RiboTaper_defined_uORFs,"Riboseq-uORFs",ifelse(TE_CDS$gene_id %in% Genes_w_only_seq_predicted_uORFs, "Predicted-uORFs", "Others"))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("Riboseq-uORFs", "Predicted-uORFs", "Others"), labels=c("Riboseq-uORFs", "Predicted-uORFs", "Others"))
TE_CDS %>% group_by(Category) %>% dplyr::summarise_at(vars(-gene_id),list(mean = mean, median = median))


TE_CDS2 <- TE_CDS %>% filter(Category %in% c("Riboseq-uORFs", "Predicted-uORFs")) %>% droplevels()
TE_CDS2$Category <- factor(TE_CDS2$Category, levels=c("Predicted-uORFs","Riboseq-uORFs"), labels=c("spuORFs","atuORF"))

library(plyr)
mu <- ddply(TE_CDS2, "Category", summarise, grp.mean=median(RNA))

pTE_CDS_hist <- ggplot(TE_CDS2,aes(x=RNA,color=Category,fill=Category)) +
  geom_histogram(alpha=0.8, position = 'identity') +
  xlim(c(0,70)) +
  theme_classic() + geom_vline(data=mu, aes(xintercept=grp.mean, color=Category),linetype = 9,size=1)
pTE_CDS_hist


pTE_CDS_fpoly <- ggplot(TE_CDS2,aes(x=RNA,color=Category)) +
  geom_freqpoly(binwidth = 1) +
  xlim(c(0,5000)) +
  theme_classic() + 
  xlab("RNA (TPM)")+
  ylab("Numbers") +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=Category),linetype = 2,size=0.5)

pTE_CDS_fpoly

# ggsave("~/Desktop/uORFs_miRNA/Figures_May052021/Figure2D_frepoly_atuORFs_vs_spuORFs.pdf",
#        plot = pTE_CDS_fpoly,
#        units = "in",
#        width = 5,
#        height = 4)

