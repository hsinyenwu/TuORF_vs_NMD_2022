#############################################################################################
# Figure 4I-L: uORF number for NMD vs Control mRNA ratio 
#############################################################################################
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

# Load NMD targets
NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
colnames(NMD) <- "gene_id"

#Load datasets
RNA <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_RNA.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
Ribo <- read.delim("~/Desktop/uORFs_miRNA/RSEM_Araport11_1st+2nd_Riboseq_Ribo.genes.results_CDS_only",header=T,sep="\t",stringsAsFactors = F, quote="")
TE <- Ribo$TPM/RNA$TPM
TE_CDS <- data.frame(gene_id=RNA$gene_id,
                     RNA=RNA$TPM,
                     Ribo=Ribo$TPM,
                     TE=TE)
# Remove not translated transcripts
TE_CDS <- TE_CDS %>% filter(gene_id %in% fiveUTRByGeneNames)
# TE_CDS <- TE_CDS %>% filter(TE>0 & !is.infinite(TE))
TE_CDS <- TE_CDS %>% filter(RNA>0)

#RiboTaper output
ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
table(ORF_max_filt$category)
# dORF     ncORFS  ORFs_ccds Overl_dORF       uORF 
# 208        379      37361         15       2093
ORF_max_filt_uORF <- ORF_max_filt %>% filter(category=="uORF")
head(ORF_max_filt_uORF,2)
nrow(ORF_max_filt_uORF)

ORF_max_filt_uORF <- ORF_max_filt_uORF %>% mutate(pept_length=nchar(ORF_pept))

atuORF_gene <- table(ORF_max_filt_uORF$gene_id)
atuORF_gene_df <- data.frame(gene_id=names(atuORF_gene),freq=as.numeric(atuORF_gene))
head(atuORF_gene_df)
atuORF_gene_df1 <- atuORF_gene_df %>% filter(freq==1)
atuORF_gene_df2 <- atuORF_gene_df %>% filter(freq>1)

#remove genes with more than one uORF from TE_CDS (so others will not include them)
# TE_CDS <- TE_CDS %>% filter(!(gene_id %in% atuORF_gene_df2$gene_id))

# ORF_max_filt_uORF <- ORF_max_filt_uORF %>% filter(gene_id %in% atuORF_gene_df1$gene_id)
# 
# #ORF_P_sites<-Total sum of P-sites position mapped to the ORF
# #ORF_RNA_sites<-Total sum of RNA-sites position mapped to the ORF
# #P_sites_sum<-Total sum of P-sites position mapped to the transcript
# #RNA_sites<-Total sum of RNA-sites position mapped to the transcript
# 
# ORF_max_filt_mORF <- ORF_max_filt %>% filter(transcript_id %in% ORF_max_filt_uORF$transcript_id,category=="ORFs_ccds")%>% group_by(transcript_id) %>% top_n(n=1, wt = ORF_length)
# length(unique(ORF_max_filt_uORF$transcript_id))
# length(unique(ORF_max_filt_mORF$transcript_id))
# ORF_max_filt_mORF <- ORF_max_filt_mORF %>% mutate(mORF_P_sites=ORF_P_sites,mORF_length=ORF_length) %>% dplyr::select(transcript_id,mORF_P_sites,mORF_length)
# ORF_max_filt_uORF_mORF <- inner_join(ORF_max_filt_uORF,ORF_max_filt_mORF) %>% mutate(uORF_mORF_Ribo_ratio=(ORF_P_sites/ORF_length)/(mORF_P_sites/mORF_length))
# dim(ORF_max_filt_uORF_mORF) #[1] 1914   42
# ORF_max_filt_uORF_low <- dplyr::top_n(ORF_max_filt_uORF_mORF,-719,uORF_mORF_Ribo_ratio)
# ORF_max_filt_uORF_high <- dplyr::top_n(ORF_max_filt_uORF_mORF,719,uORF_mORF_Ribo_ratio)

TE_CDS$Category <- ifelse(TE_CDS$gene_id %in% NMD$gene_id,"NMD",ifelse(TE_CDS$gene_id %in% atuORF_gene_df1$gene_id, "One", ifelse(TE_CDS$gene_id %in% atuORF_gene_df2$gene_id,"2orMore","Others")))
TE_CDS$Category <- factor(TE_CDS$Category, levels=c("NMD", "One", "2orMore", "Others"), labels=c("NMD", "One", "2orMore", "Others"))
table(TE_CDS$Category)
# NMD     One 2orMore  Others 
# 250    1463     240   19663 

TE_CDS %>% group_by(Category) %>% dplyr::summarise_at(vars(-gene_id),list(mean = mean, median = median))

# # A tibble: 4 × 7
# Category RNA_mean Ribo_mean TE_mean RNA_median Ribo_median TE_median
# <fct>       <dbl>     <dbl>   <dbl>      <dbl>       <dbl>     <dbl>
# 1 NMD          17.5      20.7   0.910       4.74        2.90     0.728
# 2 One          32.0      37.4   0.894      11.8         8.47     0.772
# 3 2orMore      27.6      28.7   0.671      12.8         7.32     0.583
# 4 Others       32.0      44.4   1.38        7.84        7.88     1.12 

#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#########@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

###############################
# D. Gloggnitzer pad4, upf1
###############################
# rm(list=ls())
# NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
# colnames(NMD) <- "gene_id"

pad4_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/pad4_1/abundance.tsv",header=T)
pad4_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/pad4_2/abundance.tsv",header=T)
pad4_3 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/pad4_3/abundance.tsv",header=T)
smg7pad4_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/smg7pad4_1/abundance.tsv",header=T)
smg7pad4_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/smg7pad4_2/abundance.tsv",header=T)
smg7pad4_3 <- read.delim(file="/Users/wu/Desktop/NMD paper/Gloggnitzer2014/kallisto/smg7pad4_3/abundance.tsv",header=T)
pad4_1$gene_id <- substr(pad4_1$target_id,1,9)
pad4_2$gene_id <- substr(pad4_2$target_id,1,9)
pad4_3$gene_id <- substr(pad4_3$target_id,1,9)
smg7pad4_1$gene_id <- substr(pad4_1$target_id,1,9)
smg7pad4_2$gene_id <- substr(smg7pad4_2$target_id,1,9)
smg7pad4_3$gene_id <- substr(smg7pad4_3$target_id,1,9)

pad4_1_genelevel <- pad4_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
pad4_2_genelevel <- pad4_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
pad4_3_genelevel <- pad4_3 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
smg7pad4_1_genelevel <- smg7pad4_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
smg7pad4_2_genelevel <- smg7pad4_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
smg7pad4_3_genelevel <- smg7pad4_3 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()

pad4_genelevel_df <- data.frame(gene_id=pad4_1_genelevel$gene_id,RNA=(pad4_1_genelevel$gene_tpm+pad4_1_genelevel$gene_tpm+pad4_1_genelevel$gene_tpm)/3)
smg7pad4_genelevel_df <- data.frame(gene_id=smg7pad4_1_genelevel$gene_id,RNA=(smg7pad4_1_genelevel$gene_tpm+smg7pad4_1_genelevel$gene_tpm+smg7pad4_1_genelevel$gene_tpm)/3)

#############
smg7pad4_pad4_genelevel_ratio <- inner_join(x=pad4_genelevel_df,y=smg7pad4_genelevel_df,by="gene_id")
smg7pad4_pad4_genelevel_ratio$ratio <- smg7pad4_pad4_genelevel_ratio$RNA.y/smg7pad4_pad4_genelevel_ratio$RNA.x

smg7pad4_pad4_genelevel_ratio$Category <- ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% atuORF_gene_df1$gene_id, "One", ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% atuORF_gene_df2$gene_id,"2orMore","Others")))
smg7pad4_pad4_genelevel_ratio$Category <- factor(smg7pad4_pad4_genelevel_ratio$Category, levels=c("NMD", "One", "2orMore", "Others"), labels=c("NMD", "One", "2orMore", "Others"))

# smg7pad4_pad4_genelevel_ratio$Category <- ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_low$gene_id, "Low", ifelse(smg7pad4_pad4_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_high$gene_id,"High","Others")))
# smg7pad4_pad4_genelevel_ratio$Category <- factor(smg7pad4_pad4_genelevel_ratio$Category, levels=c("NMD", "Low", "High", "Others"), labels=c("NMD", "Low", "High", "Others"))
dim(smg7pad4_pad4_genelevel_ratio) #[1] 27655     5
smg7pad4_pad4_genelevel_ratio <- smg7pad4_pad4_genelevel_ratio %>% filter(RNA.y>0 & RNA.x>0)
dim(smg7pad4_pad4_genelevel_ratio) #[1] 21144     5
smg7pad4_pad4_genelevel_ratio%>% group_by(Category) %>% dplyr::summarize(Median = median(ratio, na.rm=TRUE)/.936)
# A tibble: 4 × 2
# Category Median
# <fct>     <dbl>
#   1 NMD        2.65
# 2 One        1.08
# 3 2orMore    1.09
# 4 Others     1.00

#######################
# smg7pad4/pad4 RNA ratio boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(ratio ~ Category, data = smg7pad4_pad4_genelevel_ratio)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(ratio ~ Category, data = smg7pad4_pad4_genelevel_ratio, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- smg7pad4_pad4_genelevel_ratio %>%
  wilcox_test(ratio ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(smg7pad4_pad4_genelevel_ratio, 
                 x = "Category", y = "ratio", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.1,6.5),
                 legend = "right") + 
  geom_hline(yintercept=median(smg7pad4_pad4_genelevel_ratio$ratio[smg7pad4_pad4_genelevel_ratio$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "D") + 
  labs(title = "smg7pad4/pad4 ratio") +
  ylab("smg7pad4/pad4 RNA ratio")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(face = "italic"))

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

vsteps <- yvalue(ystart=5.2,stepsize=0.2,num=6)

psmg7pad4_pad4_RNA_ratio <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.00001, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.2,size = 2.5)
psmg7pad4_pad4_RNA_ratio

#######################################
#######################################
#######################################
# C. Raxwal 2020 upf1pad4/pad4 RNAseq
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
pad4_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/pad4_T_R1/abundance.tsv",header=T)
pad4_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/pad4_T_R2/abundance.tsv",header=T)
pad4_3 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/pad4_T_R3/abundance.tsv",header=T)
upf1pad4_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/upf1pad4_T_R1//abundance.tsv",header=T)
upf1pad4_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/upf1pad4_T_R2/abundance.tsv",header=T)
upf1pad4_3 <- read.delim(file="/Users/wu/Desktop/NMD paper/Raxwal_data/kallisto/upf1pad4_T_R3/abundance.tsv",header=T)
pad4_1$gene_id <- substr(pad4_1$target_id,1,9)
pad4_2$gene_id <- substr(pad4_2$target_id,1,9)
pad4_3$gene_id <- substr(pad4_3$target_id,1,9)
upf1pad4_1$gene_id <- substr(pad4_1$target_id,1,9)
upf1pad4_2$gene_id <- substr(upf1pad4_2$target_id,1,9)
upf1pad4_3$gene_id <- substr(upf1pad4_3$target_id,1,9)

pad4_1_genelevel <- pad4_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
pad4_2_genelevel <- pad4_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
pad4_3_genelevel <- pad4_3 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf1pad4_1_genelevel <- upf1pad4_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf1pad4_2_genelevel <- upf1pad4_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf1pad4_3_genelevel <- upf1pad4_3 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()

pad4_genelevel_df <- data.frame(gene_id=pad4_1_genelevel$gene_id,RNA=(pad4_1_genelevel$gene_tpm+pad4_1_genelevel$gene_tpm+pad4_1_genelevel$gene_tpm)/3) %>% filter(RNA>0)
upf1pad4_genelevel_df <- data.frame(gene_id=upf1pad4_1_genelevel$gene_id,RNA=(upf1pad4_1_genelevel$gene_tpm+upf1pad4_1_genelevel$gene_tpm+upf1pad4_1_genelevel$gene_tpm)/3) %>% filter(RNA>0)

upf1pad4_pad4_genelevel_ratio <- inner_join(x=pad4_genelevel_df,y=upf1pad4_genelevel_df,by="gene_id")
upf1pad4_pad4_genelevel_ratio$ratio <- upf1pad4_pad4_genelevel_ratio$RNA.y/upf1pad4_pad4_genelevel_ratio$RNA.x

upf1pad4_pad4_genelevel_ratio$Category <- ifelse(upf1pad4_pad4_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(upf1pad4_pad4_genelevel_ratio$gene_id %in% atuORF_gene_df1$gene_id, "One", ifelse(upf1pad4_pad4_genelevel_ratio$gene_id %in% atuORF_gene_df2$gene_id,"2orMore","Others")))
upf1pad4_pad4_genelevel_ratio$Category <- factor(upf1pad4_pad4_genelevel_ratio$Category, levels=c("NMD", "One", "2orMore", "Others"), labels=c("NMD", "One", "2orMore", "Others"))

dim(upf1pad4_pad4_genelevel_ratio) #[1] 21828     5
upf1pad4_pad4_genelevel_ratio <- upf1pad4_pad4_genelevel_ratio %>% filter(RNA.y>0 & RNA.x>0)
dim(upf1pad4_pad4_genelevel_ratio) #[1] 21828     5

table(upf1pad4_pad4_genelevel_ratio$Category)
upf1pad4_pad4_genelevel_ratio%>% group_by(Category) %>% dplyr::summarize(Median = median(ratio, na.rm=TRUE)/0.872)

# # A tibble: 4 × 2
# Category Median
# <fct>     <dbl>
#   1 NMD       4.99 
# 2 One       1.07 
# 3 2orMore   0.993
# 4 Others    1.00 

#######################
# upf1pad4/pad4 RNA ratio boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(ratio ~ Category, data = upf1pad4_pad4_genelevel_ratio)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(ratio ~ Category, data = upf1pad4_pad4_genelevel_ratio, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- upf1pad4_pad4_genelevel_ratio %>%
  wilcox_test(ratio ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(upf1pad4_pad4_genelevel_ratio, 
                 x = "Category", y = "ratio", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.3,17),
                 legend = "right") + 
  geom_hline(yintercept=median(upf1pad4_pad4_genelevel_ratio$ratio[upf1pad4_pad4_genelevel_ratio$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "C") + 
  labs(title = "upf1pad4/pad4 ratio") +
  ylab("upf1pad4/pad4 RNA ratio")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(face = "italic"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

vsteps <- yvalue(ystart=13,stepsize=0.6,num=6)

pupf1pad4_pad4_RNA_ratio <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0000000005, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.2,size = 2.5)
pupf1pad4_pad4_RNA_ratio
# psmg7pad4_pad4_RNA_ratio



#######################
# B. Merchante 2015 upf2/Col RNAseq
#######################
# rm(list=ls())
ColRNA1 <- read.delim(file="~/Desktop/Merchante2015/Kallisto/kallisto/Col-Air-Total1/abundance.tsv",header=T)
ColRNA2 <- read.delim(file="~/Desktop/Merchante2015/Kallisto/kallisto/Col-Air-Total1/abundance.tsv",header=T)

upf2RNA1 <- read.delim(file="~/Desktop/Merchante2015/Kallisto/kallisto/upf2-Air-Total1/abundance.tsv",header=T)
upf2RNA2 <- read.delim(file="~/Desktop/Merchante2015/Kallisto/kallisto/upf2-Air-Total1/abundance.tsv",header=T)

ColRNA1$gene_id <- substr(ColRNA1$target_id,1,9)
ColRNA2$gene_id <- substr(ColRNA2$target_id,1,9)
upf2RNA1$gene_id <- substr(upf2RNA1$target_id,1,9)
upf2RNA2$gene_id <- substr(upf2RNA2$target_id,1,9)

ColRNA1_genelevel <- ColRNA1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
ColRNA2_genelevel <- ColRNA2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf2RNA1_genelevel <- upf2RNA1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf2RNA2_genelevel <- upf2RNA2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()

Col_genelevel_df <- inner_join(ColRNA1_genelevel,ColRNA2_genelevel,by="gene_id")%>% mutate(RNA=(gene_tpm.x+gene_tpm.y)/2)%>% filter(RNA>0 ) %>% filter(gene_id %in% fiveUTRByGeneNames)
upf2_genelevel_df <- inner_join(upf2RNA1_genelevel,upf2RNA2_genelevel,by="gene_id")%>% mutate(RNA=(gene_tpm.x+gene_tpm.y)/2)%>% filter(RNA>0 ) %>% filter(gene_id %in% fiveUTRByGeneNames) 

##############
# Col_genelevel_df$Category <- ifelse(Col_genelevel_df$gene_id %in% MuORF_no_TuORF_gid,"MuORFs",ifelse(Col_genelevel_df$gene_id %in% TuORF_no_MuORF_gid, "TuORFs", ifelse(Col_genelevel_df$gene_id %in% MuORF_and_TuORF_gid,"Both","Others")))
# Col_genelevel_df$Category <- factor(Col_genelevel_df$Category, levels=c("MuORFs", "TuORFs", "Both", "Others"), labels=c("MuORFs", "TuORFs", "Both", "Others"))
# 
# upf2_genelevel_df$Category <- ifelse(upf2_genelevel_df$gene_id %in% MuORF_no_TuORF_gid,"MuORFs",ifelse(upf2_genelevel_df$gene_id %in% TuORF_no_MuORF_gid, "TuORFs", ifelse(upf2_genelevel_df$gene_id %in% MuORF_and_TuORF_gid,"Both","Others")))
# upf2_genelevel_df$Category <- factor(upf2_genelevel_df$Category, levels=c("MuORFs", "TuORFs", "Both", "Others"), labels=c("MuORFs", "TuORFs", "Both", "Others"))
# 
upf2_Col_genelevel_ratio <- inner_join(x=Col_genelevel_df,y=upf2_genelevel_df,by="gene_id") %>% filter(gene_id %in% fiveUTRByGeneNames)
upf2_Col_genelevel_ratio$ratio <- upf2_Col_genelevel_ratio$RNA.y/upf2_Col_genelevel_ratio$RNA.x

upf2_Col_genelevel_ratio$Category <- ifelse(upf2_Col_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(upf2_Col_genelevel_ratio$gene_id %in% atuORF_gene_df1$gene_id, "One", ifelse(upf2_Col_genelevel_ratio$gene_id %in% atuORF_gene_df2$gene_id,"2orMore","Others")))
upf2_Col_genelevel_ratio$Category <- factor(upf2_Col_genelevel_ratio$Category, levels=c("NMD", "One", "2orMore", "Others"), labels=c("NMD", "One", "2orMore", "Others"))


# upf2_Col_genelevel_ratio$Category <- ifelse(upf2_Col_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(upf2_Col_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_low$gene_id, "Low", ifelse(upf2_Col_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_high$gene_id,"High","Others")))
# upf2_Col_genelevel_ratio$Category <- factor(upf2_Col_genelevel_ratio$Category, levels=c("NMD", "Low", "High", "Others"), labels=c("NMD", "Low", "High", "Others"))
dim(upf2_Col_genelevel_ratio) #[1] 20722     5
upf2_Col_genelevel_ratio <- upf2_Col_genelevel_ratio %>% filter(RNA.y>0 & RNA.x>0)
dim(upf2_Col_genelevel_ratio) #[1] 20722     

# upf2_Col_genelevel_ratio$Category <- ifelse(upf2_Col_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_low$gene_id,"Low",ifelse(upf2_Col_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_mid$gene_id, "Mid", ifelse(upf2_Col_genelevel_ratio$gene_id %in% ORF_max_filt_uORF_high$gene_id,"High","Others")))
# upf2_Col_genelevel_ratio$Category <- factor(upf2_Col_genelevel_ratio$Category, levels=c("Low", "Mid", "High", "Others"), labels=c("Low", "Mid", "High", "Others"))

# upf2_Col_genelevel_ratio$Category <- ifelse(upf2_Col_genelevel_ratio$gene_id %in% NMD$gene_id, "NMD targets", ifelse(upf2_Col_genelevel_ratio$gene_id %in% MuORF_no_TuORF_gid,"MuORFs",ifelse(upf2_Col_genelevel_ratio$gene_id %in% TuORF_no_MuORF_gid, "TuORFs", ifelse(upf2_Col_genelevel_ratio$gene_id %in% MuORF_and_TuORF_gid,"MuORFs&TuORFs","Others"))))
# upf2_Col_genelevel_ratio$Category <- factor(upf2_Col_genelevel_ratio$Category, levels=c("NMD targets", "MuORFs", "TuORFs", "MuORFs&TuORFs", "Others"), labels=c("NMD targets", "MuORFs", "TuORFs", "MuORFs&TuORFs", "Others"))
# upf2_Col_genelevel_ratio <- upf2_Col_genelevel_ratio %>% filter(!(gene_id %in% NMD$gene_id))


upf2_Col_genelevel_ratio %>% group_by(Category) %>% dplyr::summarize(Median = median(ratio, na.rm=TRUE)/0.888)
# A tibble: 4 × 2
# Category Median
# <fct>     <dbl>
# 1 NMD       1.72 
# 2 One       0.983
# 3 2orMore   0.976
# 4 Others    1.00
#######################
# upf2/Col RNA ratio boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(ratio ~ Category, data = upf2_Col_genelevel_ratio)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(ratio ~ Category, data = upf2_Col_genelevel_ratio, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- upf2_Col_genelevel_ratio %>%
  wilcox_test(ratio ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(upf2_Col_genelevel_ratio, 
                 x = "Category", y = "ratio", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.5,6),
                 legend = "right") + 
  geom_hline(yintercept=median(upf2_Col_genelevel_ratio$ratio[upf2_Col_genelevel_ratio$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "B") + 
  labs(title = "upf2/Col ratio") +
  ylab("upf2/Col RNA ratio")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(face = "italic"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

vsteps <- yvalue(ystart=4.5,stepsize=0.26,num=6)
pupf2_Col_RNA_ratio <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.2,size = 2.5)
pupf2_Col_RNA_ratio

#######################
#######################
#######################
#######################
#######################
# A. upf1upf3/Col RNAseq
#######################

Col_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Drechsel2013/kallisto/NMD_WT_R1/abundance.tsv",header=T)
Col_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Drechsel2013/kallisto/NMD_WT_R2/abundance.tsv",header=T)

upf1upf3_1 <- read.delim(file="/Users/wu/Desktop/NMD paper/Drechsel2013/kallisto/NMD_UPF1UPF3_R1/abundance.tsv",header=T)
upf1upf3_2 <- read.delim(file="/Users/wu/Desktop/NMD paper/Drechsel2013/kallisto/NMD_UPF1UPF3_R2/abundance.tsv",header=T)

Col_1$gene_id <- substr(Col_1$target_id,1,9)
Col_2$gene_id <- substr(Col_2$target_id,1,9)

upf1upf3_1$gene_id <- substr(Col_1$target_id,1,9)
upf1upf3_2$gene_id <- substr(upf1upf3_2$target_id,1,9)

Col_1_genelevel <- Col_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
Col_2_genelevel <- Col_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf1upf3_1_genelevel <- upf1upf3_1 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()
upf1upf3_2_genelevel <- upf1upf3_2 %>% group_by(gene_id) %>% summarise(gene_tpm = sum(tpm)) %>% as.data.frame()

Col_genelevel_df <- data.frame(gene_id=Col_1_genelevel$gene_id,RNA=(Col_1_genelevel$gene_tpm+Col_2_genelevel$gene_tpm)/2) %>% filter(RNA>0)
upf1upf3_genelevel_df <- data.frame(gene_id=upf1upf3_1_genelevel$gene_id,RNA=(upf1upf3_1_genelevel$gene_tpm+upf1upf3_2_genelevel$gene_tpm)/2) %>% filter(RNA>0)

upf1upf3_Col_genelevel_ratio <- inner_join(x=Col_genelevel_df,y=upf1upf3_genelevel_df,by="gene_id")
upf1upf3_Col_genelevel_ratio$ratio <- upf1upf3_Col_genelevel_ratio$RNA.y/upf1upf3_Col_genelevel_ratio$RNA.x
upf1upf3_Col_genelevel_ratio$Category <- ifelse(upf1upf3_Col_genelevel_ratio$gene_id %in% NMD$gene_id,"NMD",ifelse(upf1upf3_Col_genelevel_ratio$gene_id %in% atuORF_gene_df1$gene_id, "One", ifelse(upf1upf3_Col_genelevel_ratio$gene_id %in% atuORF_gene_df2$gene_id,"2orMore","Others")))
upf1upf3_Col_genelevel_ratio$Category <- factor(upf1upf3_Col_genelevel_ratio$Category, levels=c("NMD", "One", "2orMore", "Others"), labels=c("NMD", "One", "2orMore", "Others"))
dim(upf1upf3_Col_genelevel_ratio) #[1] 20722     5
upf1upf3_Col_genelevel_ratio <- upf1upf3_Col_genelevel_ratio %>% filter(RNA.y>0 & RNA.x>0)
dim(upf1upf3_Col_genelevel_ratio) #[1] 20722     

upf1upf3_Col_genelevel_ratio %>% group_by(Category) %>% dplyr::summarize(Median = median(ratio, na.rm=TRUE)/1.06)

# # A tibble: 4 × 2
# Category Median
# <fct>     <dbl>
# 1 NMD       3.05 
# 2 One       1.02 
# 3 2orMore   0.929
# 4 Others    0.998
#######################
# upf1upf3/Col RNA ratio boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(ratio ~ Category, data = upf1upf3_Col_genelevel_ratio)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(ratio ~ Category, data = upf1upf3_Col_genelevel_ratio, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- upf1upf3_Col_genelevel_ratio %>%
  wilcox_test(ratio ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(upf1upf3_Col_genelevel_ratio, 
                 x = "Category", y = "ratio", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.3,15),
                 legend = "right") + 
  geom_hline(yintercept=median(upf1upf3_Col_genelevel_ratio$ratio[upf1upf3_Col_genelevel_ratio$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "A") + 
  labs(title = "upf1upf3/Col ratio") +
  ylab("RNA (TPM)")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(face = "italic"))

# Add p-values onto the box plots
stat.test <- stat.test %>%
  add_xy_position(x = "Category", dodge = 0.8)

vsteps <- yvalue(ystart=12.35,stepsize=0.5,num=6)

pupf1upf3_Col_RNA_ratio <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.0000000003, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.2,size = 2.5)
pupf1upf3_Col_RNA_ratio

##########################################
##########################################
##########################################
get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pupf2_Col_RNA_ratio)

# pupf1upf3_Col_RNA_ratio
# pupf2_Col_RNA_ratio
# pupf1pad4_pad4_RNA_ratio
# psmg7pad4_pad4_RNA_ratio

p4 <- grid.arrange(pupf1upf3_Col_RNA_ratio + theme(legend.position="none"),
                   pupf2_Col_RNA_ratio + theme(legend.position="none"),
                   pupf1pad4_pad4_RNA_ratio + theme(legend.position="none"),
                   psmg7pad4_pad4_RNA_ratio + theme(legend.position="none"),
                   p2_legend,
                   layout_matrix=rbind(c(1,2,3,4,5)),
                   top=textGrob("   Figure x. TuORF Ribo compared to main ORF Ribo ratio vs mRNA ratio betweenNMD mutants and WT", gp=gpar(fontsize=13), x = 0, hjust = 0))

ggsave("~/Desktop/atRTD3/TuORF_number_vs_NMDt_vs_Other_in_NMD_mutants.pdf",
       plot = p4,
       units = "in",
       width = 14.5,
       height = 5)
```
