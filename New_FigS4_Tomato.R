################
# Figure 2E-F Tomato
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


# ggsave("~/Desktop/FigureS4_tomato_RNA_TE_box.pdf",
#        plot = p3,
#        units = "in",
#        width = 7,
#        height = 3)
