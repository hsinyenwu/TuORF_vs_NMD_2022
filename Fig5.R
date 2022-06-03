########################################################################
# Figure 4. NMD targets, uORF genes, other genes RNA TE Half-lives Protein
########################################################################
#
rm(list=ls())
library(openxlsx)
library(dplyr)
library(reshape)
library(ggplot2)
library(ggtext)
library(grid)
library(gridExtra)
library(ggpubr)
library(rstatix)
library(EnvStats)

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
library(clipr)
write_clip(TE_CDS_NMD$gene_id)

# Expression levels (for plotting RNA-seq and TE)
TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(TE_CDS$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(TE_CDS$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))
table(TE_CDS$Category)

# mRNA half-lives
half_lives <- read.delim(file="/Volumes/GoogleDrive/My Drive/CPuORF_paper/data-for-analysis/RNA-half-lives-5EU.txt",header=T,stringsAsFactors=F,sep="\t")
half_lives <- half_lives %>% dplyr::select("Gene.ID","Mean")
colnames(half_lives) <- c("gene_id","Mean")
half_lives <- half_lives %>% filter(gene_id %in% fiveUTRByGeneNames)
half_lives$Category <- ifelse(half_lives$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(half_lives$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(half_lives$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
half_lives$Category <- factor(half_lives$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))

# Quantitative peoteomics for Shoot
shoot_proteomics <- read.xlsx("~/Desktop/uORFs_miRNA/Justin_Proteomics2017/data_excel/pmic12946-sup-0004-tables7.xlsx",sheet = 3)
head(shoot_proteomics$`MS/MS.count`)
shoot_proteomics$gene_id <- substr(shoot_proteomics$Protein.IDs,1,9)
shoot_proteomics2 <- shoot_proteomics %>% dplyr::select(gene_id,`MS/MS.count`,`Sequence.length`)
colnames(shoot_proteomics2) <-c("gene_id","MS_count","Pept_length")
shoot_proteomics2 <- shoot_proteomics2 %>% filter(gene_id %in% fiveUTRByGeneNames)
shoot_proteomics2$norm_MS_count <- shoot_proteomics2$MS_count/shoot_proteomics2$Pept_length
shoot_proteomics2$Category <- ifelse(shoot_proteomics2$gene_id %in% NMD_target_wo_uORFs,"NMD_target_wo_uORFs",ifelse(shoot_proteomics2$gene_id %in% NMD_target_w_uORFs, "NMD_target_w_uORFs", ifelse(shoot_proteomics2$gene_id %in% uORFs_not_NMD_targets,"uORFs_not_NMD_targets","Others")))
shoot_proteomics2$Category <- factor(shoot_proteomics2$Category, levels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"), labels=c("NMD_target_wo_uORFs", "NMD_target_w_uORFs", "uORFs_not_NMD_targets", "Others"))

#
TE_CDS %>% group_by(Category) %>% summarise(median=median(RNA)/8.19)
TE_CDS %>% group_by(Category) %>% summarise(median=median(TE)/1.14)
half_lives %>% group_by(Category) %>% summarise(median=median(Mean)/1.87)
shoot_proteomics2 %>% group_by(Category) %>% summarise(median=median(norm_MS_count)/0.0440)

pRNA <- ggplot(TE_CDS, aes(x=RNA,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,50)+xlab("RNA (TPM)")+ylab("Fraction of Genes")+labs(tag = "A") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12))
pRNA 

pTE <- ggplot(TE_CDS, aes(x=TE,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,4)+xlab("TE")+ylab("Fraction of Genes")+labs(tag = "B") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12))
pTE 


pHL <- ggplot(half_lives, aes(x=Mean,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,5)+xlab("Half-lives (Hours)")+ylab("Fraction of Genes")+labs(tag = "C") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12))
pHL 

pProtein <- ggplot(shoot_proteomics2, aes(x=norm_MS_count,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.5)+xlab("Half-lives (Hours)")+ylab("Fraction of Genes")+labs(tag = "D") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12))
pProtein

#a function to get legend from p_mdecay_uORFs
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pProtein)

p3= grid.arrange(pRNA + theme(legend.position="none"),
                 pTE + theme(legend.position="none"),
                 pHL + theme(legend.position="none"),
                 pProtein + theme(legend.position="none"),
                 p2_legend,
                 layout_matrix=rbind(c(1,2,3,4,5)),
                 top=textGrob("Figure 5. ECDF",gp = gpar(fontsize = 12))) #fontface = "bold"

p3

# dir.create("~/Desktop/NMD_revise/Figures/temp")
# ggsave("~/Desktop/NMD_revise/Figures/temp/Figure5ECDF.pdf",
#        plot = p3,
#        units = "in",
#        width = 13,
#        height = 2.7)


#######################
####   Box plot    ####
#######################

#######################
# RNA boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(RNA ~ Category, data = TE_CDS)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(RNA ~ Category, data = TE_CDS, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- TE_CDS %>%
  wilcox_test(RNA ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(TE_CDS, 
                 x = "Category", y = "RNA", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-1,80),
                 legend = "right") + 
  geom_hline(yintercept=median(TE_CDS$RNA[TE_CDS$Category=="uORFs_not_NMD_targets"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "A") + 
  ylab("RNA (TPM)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

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

vsteps <- yvalue(ystart=53,stepsize=4,num=6)

pRNA_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.00002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-1.5,size = 2.5)
pRNA_box

#######################
# TE boxplot
#######################

TE_CDS <- TE_CDS[!is.infinite(TE_CDS$TE),]
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(TE ~ Category, data = TE_CDS)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(TE ~ Category, data = TE_CDS, var.equal = FALSE)
# p-value = 2.2e-16

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
                 ylim=c(-0.1,4.1),
                 legend = "right") + 
  geom_hline(yintercept=median(TE_CDS$TE[TE_CDS$Category=="uORFs_not_NMD_targets"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "C") + 
  ylab("TE (TPM)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12))

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

vsteps <- yvalue(ystart=2.8,stepsize=0.22,num=6)

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

pTE_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.00002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.1,size = 2.5)
pTE_box

#######################
# Half-life boxplot
#######################
half_lives <- half_lives[!is.infinite(half_lives$Mean),]
#test equal variance
bartlett.test(Mean ~ Category, data = half_lives)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(Mean ~ Category, data = half_lives, var.equal = FALSE)
# p-value = 2.2e-16

# default p-value adjustement method: Holm
stat.test <- half_lives %>%
  t_test(Mean ~Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

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
                 notch = F,
                 ylim=c(-0.1,6),
                 legend = "right",
                 ylab="Mean half-life") + 
  geom_hline(yintercept=median(half_lives$Mean[half_lives$Category=="uORFs_not_NMD_targets"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "D") + 
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) #legend key font size

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.05)

yvalue <- function(ystart,stepsize,num) {
  Ev <- c()
  for(i in 1:num){
    Ev[i] <- ystart + stepsize*(i-1)
  }
  return(Ev)
}

vsteps <- yvalue(ystart=4,stepsize=0.35,num=6)

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
                 notch = F,
                 ylim=c(-0.01,0.38),
                 legend = "right",
                 ylab="Normalized MS count") + 
  geom_hline(yintercept=median(shoot_proteomics2$norm_MS_count[shoot_proteomics2$Category=="uORFs_not_NMD_targets"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "E") + 
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) #legend key font size


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

vsteps <- yvalue(ystart=0.25,stepsize=0.026,num=6)

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

p3 = grid.arrange(pRNA_box + theme(legend.position="none"),
                  pTE_box + theme(legend.position="none"),
                  pHL_box + theme(legend.position="none"),
                  pProtein_box+ theme(legend.position="none"),
                  p2_legend,
                  layout_matrix=rbind(c(1,2,3,4,5)),
                  top=textGrob("Figure 4. NMD targets vs uORF-containing genes vs Other mRNAs: RNA/TE/Half-lives/Protein", gp=gpar(fontsize=12), x = 0, hjust = 0))

ggsave("~/Desktop/uORFs_miRNA/Figures_Aug252021/Figure 4. NMD targets vs uORF-containing genes vs Other mRNAs expression TE half-lives protein.pdf",
       plot = p3,
       units = "in",
       width = 13,
       height = 4)
       
