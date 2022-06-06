#############################################
# Z: Decay rates in different mutants vs NMD targets, uORF genes, other genes
#############################################
rm(list=ls())
library(dbplyr)
library(reshape)
library(ggplot2)
library(ggtext)
library(openxlsx)

#Load NMD gene list, mRNA decay rates, RiboTaper output
NMD <- read.xlsx("~/Desktop/uORFs_miRNA/NMD_targets_TPC2020.xlsx",sheet=5)
colnames(NMD) <- "gene_id"

decay <- read.xlsx("~/Desktop/uORFs_miRNA/Sorensona_et_al_pnas2017.xlsx",sheet=3)

# ORF_max_filt <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
Exp <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
Exp_uORFs <- Exp %>% filter(category=="uORF")
Exp_uORFs$pept_length <- nchar(Exp_uORFs$ORF_pept)
Exp_uORFs$dist_to_main <- Exp_uORFs$annotated_start-Exp_uORFs$stop_pos
Exp_uORFs_NMD <- Exp_uORFs %>% filter(gene_id %in% NMD$gene_id)

NMD <- NMD %>% filter(!(gene_id %in% Exp_uORFs$gene_id))
# Exp_uORFs_NMD
Exp_uORFs_not_NMD <- Exp_uORFs %>% filter(!(gene_id %in% NMD$gene_id))
wilcox.test(x=Exp_uORFs_not_NMD$pept_length, y=Exp_uORFs_NMD$pept_length)
t.test(x=Exp_uORFs_not_NMD$pept_length, y=Exp_uORFs_NMD$pept_length,var.equal = T)
t.test(x=Exp_uORFs_not_NMD$pept_length, y=Exp_uORFs_NMD$pept_length,var.equal = F)
median(Exp_uORFs_not_NMD$pept_length) #20.5
median(Exp_uORFs_NMD$pept_length) #28
# table(Exp_uORFs_not_NMD$gene_id)
# table(Exp_uORFs_NMD$gene_id)

median(Exp_uORFs_not_NMD$dist_to_main) #99
median(Exp_uORFs_NMD$dist_to_main) #128
mean(Exp_uORFs_not_NMD$dist_to_main) #166.2796
mean(Exp_uORFs_NMD$dist_to_main) #196.1692

median(Exp_uORFs_not_NMD$n_exons) #6
median(Exp_uORFs_NMD$n_exons) #5

#melt decay rate 
mdecay <- melt(decay, id="gene_id")
#organizing uORF dataset
WT_uORF <- mdecay %>% filter(variable=="WT") %>% filter(gene_id %in% Exp_uORFs$gene_id)
sov_uORF <- mdecay %>% filter(variable=="sov") %>% filter(gene_id %in% Exp_uORFs$gene_id)
vcs_uORF <- mdecay %>% filter(variable=="vcs") %>% filter(gene_id %in% Exp_uORFs$gene_id)
vcs.sov_uORF <- mdecay %>% filter(variable=="vcs.sov") %>% filter(gene_id %in% Exp_uORFs$gene_id)
WT_uORF$Category <- "uORF-containing genes"
sov_uORF$Category <- "uORF-containing genes"
vcs_uORF$Category <- "uORF-containing genes"
vcs.sov_uORF$Category <- "uORF-containing genes"

#organizing NMD dataset
WT_NMD <- mdecay %>% filter(variable=="WT") %>% filter(gene_id %in% NMD$gene_id)
sov_NMD <- mdecay %>% filter(variable=="sov") %>% filter(gene_id %in% NMD$gene_id)
vcs_NMD <- mdecay %>% filter(variable=="vcs") %>% filter(gene_id %in% NMD$gene_id)
vcs.sov_NMD <- mdecay %>% filter(variable=="vcs.sov") %>% filter(gene_id %in% NMD$gene_id)
WT_NMD$Category <- "NMD targets"
sov_NMD$Category <- "NMD targets"
vcs_NMD$Category <- "NMD targets"
vcs.sov_NMD$Category <- "NMD targets"

#organizing NMD and uORF dataset
WT_NMDt_TuORF <- mdecay %>% filter(variable=="WT") %>% filter(gene_id %in% Exp_uORFs_NMD$gene_id)
sov_NMDt_TuORF <- mdecay %>% filter(variable=="sov") %>% filter(gene_id %in% Exp_uORFs_NMD$gene_id)
vcs_NMDt_TuORF <- mdecay %>% filter(variable=="vcs") %>% filter(gene_id %in% Exp_uORFs_NMD$gene_id)
vcs.sov_NMDt_TuORF <- mdecay %>% filter(variable=="vcs.sov") %>% filter(gene_id %in% Exp_uORFs_NMD$gene_id)
WT_NMDt_TuORF$Category <- "NMDt with TuORFs"
sov_NMDt_TuORF$Category <- "NMDt with TuORFs"
vcs_NMDt_TuORF$Category <- "NMDt with TuORFs"
vcs.sov_NMDt_TuORF$Category <- "NMDt with TuORFs"

#organizing other genes dataset
WT_Others <- mdecay %>% filter(variable=="WT") %>% filter(!(gene_id %in% Exp_uORFs$gene_id) & !(gene_id %in% NMD$gene_id))
sov_Others <- mdecay %>% filter(variable=="sov") %>% filter(!(gene_id %in% Exp_uORFs$gene_id) & !(gene_id %in% NMD$gene_id))
vcs_Others <- mdecay %>% filter(variable=="vcs") %>% filter(!(gene_id %in% Exp_uORFs$gene_id) & !(gene_id %in% NMD$gene_id))
vcs.sov_Others <- mdecay %>% filter(variable=="vcs.sov") %>% filter(!(gene_id %in% Exp_uORFs$gene_id) & !(gene_id %in% NMD$gene_id))
WT_Others$Category <- "Others"
sov_Others$Category <- "Others"
vcs_Others$Category <- "Others"
vcs.sov_Others$Category <- "Others"

WT_df <- rbind(WT_Others,WT_uORF,WT_NMD,WT_NMDt_TuORF)
sov_df <- rbind(sov_Others,sov_uORF,sov_NMD,sov_NMDt_TuORF)
vcs_df <- rbind(vcs_Others,vcs_uORF,vcs_NMD,vcs_NMDt_TuORF)
vcs.sov_df <- rbind(vcs.sov_Others,vcs.sov_uORF,vcs.sov_NMD,vcs.sov_NMDt_TuORF)

WT_df$Category <- factor(WT_df$Category, levels=c("NMD targets", "uORF-containing genes", "NMDt with TuORFs", "Others"), labels=c("NMD targets", "uORF genes", "NMDt_TuORF", "Others"))
sov_df$Category <- factor(sov_df$Category, levels=c("NMD targets", "uORF-containing genes", "NMDt with TuORFs", "Others"), labels=c("NMD targets", "uORF genes", "NMDt_TuORF", "Others"))
vcs_df$Category <- factor(vcs_df$Category, levels=c("NMD targets", "uORF-containing genes", "NMDt with TuORFs", "Others"), labels=c("NMD targets", "uORF genes", "NMDt_TuORF", "Others"))
vcs.sov_df$Category <- factor(vcs.sov_df$Category, levels=c("NMD targets", "uORF-containing genes", "NMDt with TuORFs", "Others"), labels=c("NMD targets", "uORF genes", "NMDt_TuORF", "Others"))

WT_df %>% group_by(Category) %>% dplyr::summarise(median=median(value)/0.00612)
sov_df %>% group_by(Category) %>% dplyr::summarise(median=median(value)/0.00675)
vcs_df %>% group_by(Category) %>% dplyr::summarise(median=median(value)/0.00481)
vcs.sov_df %>% group_by(Category) %>% dplyr::summarise(median=median(value)/0.00394)

pWT <- ggplot(WT_df, aes(x=value,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.03)+xlab("Decay Rate")+ylab("Fraction of Genes")+labs(tag = "A") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12),
                          legend.text = element_text(size=12)) + 
  ggtitle("Col-0") +
  scale_colour_discrete("Category",labels=c("NMD targets","uORF genes","NMDt with TuORFs","Others"))
pWT 

psov <- ggplot(sov_df, aes(x=value,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.03)+xlab("Decay Rate")+ylab("Fraction of Genes")+labs(tag = "B") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12)) + 
  ggtitle("sov") +
  scale_colour_discrete("Category",labels=c("NMD targets","uORF genes","NMDt with TuORFs","Others"))
psov

pvcs <- ggplot(vcs_df, aes(x=value,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.03)+xlab("Decay Rate")+ylab("Fraction of Genes")+labs(tag = "C") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12)) + 
  ggtitle("vcs") +
  scale_colour_discrete("Category",labels=c("NMD targets","uORF genes","NMDt with TuORFs","Others"))
pvcs

pvcs.sov <- ggplot(vcs.sov_df, aes(x=value,color=Category))+
  stat_ecdf(geom = "step")+
  xlim(0,0.03)+xlab("Decay Rate")+ylab("Fraction of Genes")+labs(tag = "D") +
  theme_classic() + theme(text = element_text(size=12), 
                          plot.title = element_text(hjust = 0.5, size = 12, face = "italic"),
                          legend.text = element_text(size=12),
                          legend.title = element_text(size=12)) + 
  ggtitle("vcs.sov") +
  scale_colour_discrete("Category",labels=c("NMD targets","uORF genes","NMDt with TuORFs","Others"))
pvcs.sov

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pvcs.sov)

p3 =grid.arrange(pWT + theme(legend.position="none"),
                 psov + theme(legend.position="none"),
                 pvcs + theme(legend.position="none"),
                 pvcs.sov + theme(legend.position="none"),
                 p2_legend,
                 layout_matrix=rbind(c(1,2,3,4,5)),
                 top=textGrob("mRNA decay rates for NMD targets, uORF genes and others in WT and mutants",gp = gpar(fontsize = 14, fontface = "bold")))


# 
ggsave("~/Desktop/NMD_revise/Figures/Figure 6X. mRNA decay rates for uORF genes NMD-targets in wt and decay mutants_NMDt_rm_TuORFs.pdf",
       plot = p3,
       units = "in",
       width = 13,
       height = 3.5)


#######################
####   Box plot    ####
#######################

#######################
# WT boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(value ~ Category, data = WT_df)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(value ~ Category, data = WT_df, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- WT_df %>%
  wilcox_test(value ~ Category) %>% 
  adjust_pvalue(method = "BY")%>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(WT_df, 
                 x = "Category", y = "value", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.002,0.040),
                 legend = "right") + 
  geom_hline(yintercept=median(WT_df$value[WT_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "A",title= "WT") + 
  ylab("Decay rate")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 12))

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

vsteps <- yvalue(ystart=0.028,stepsize=0.0017,num=6)

pWT_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.001,size = 2.5)
pWT_box

#######################
# sov boxplot
#######################

#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(value ~ Category, data = sov_df)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(value ~ Category, data = sov_df, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- sov_df %>%
  wilcox_test(value ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(sov_df, 
                 x = "Category", y = "value", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.002,0.040),
                 legend = "right") + 
  geom_hline(yintercept=median(sov_df$value[sov_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "B",title= "sov") + 
  ylab("Decay rate")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))

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

vsteps <- yvalue(ystart=0.028,stepsize=0.0017,num=6)

psov_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.001,size = 2.5)
psov_box


#######################
# vcs boxplot
#######################

#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(value ~ Category, data = vcs_df)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(value ~ Category, data = vcs_df, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- vcs_df %>%
  wilcox_test(value ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(vcs_df, 
                 x = "Category", y = "value", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.002,0.040),
                 legend = "right") + 
  geom_hline(yintercept=median(vcs_df$value[vcs_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "C",title= "vcs") + 
  ylab("Decay rate")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))

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

vsteps <- yvalue(ystart=0.028,stepsize=0.0017,num=6)

pvcs_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.001,size = 2.5)
pvcs_box

#######################
# vcs.sov boxplot
#######################
#http://www.biostathandbook.com/kruskalwallis.html
#test equal variance
bartlett.test(value ~ Category, data = vcs.sov_df)
# p-value < 2.2e-16 -> reject the hypothesis for equal variance
# Welch's anova for data with unequal variance
oneway.test(value ~ Category, data = vcs.sov_df, var.equal = FALSE)
# p-value = 2.2e-16

stat.test <- vcs.sov_df %>%
  wilcox_test(value ~ Category) %>% 
  adjust_pvalue(method = "BY") %>%
  add_significance("p.adj")
stat.test

# Create a box plot
bxp <- ggboxplot(vcs.sov_df, 
                 x = "Category", y = "value", 
                 color = "Category", 
                 outlier.shape=NA,
                 notch = FALSE,
                 ylim=c(-0.002,0.040),
                 legend = "right") + 
  geom_hline(yintercept=median(vcs.sov_df$value[vcs.sov_df$Category=="Others"],na.rm=T), linetype="dashed", color = "grey30", size=0.5) +
  rotate_x_text(angle = 45) + 
  labs(tag = "D",title= "vcs.sov") + 
  ylab("Decay rate")+
  theme(plot.tag = element_text(face = 'bold',size=12),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size=12),
        axis.text.x=element_blank(),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5, size = 12, face = "italic"))

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

vsteps <- yvalue(ystart=0.028,stepsize=0.0017,num=6)

pvcs.sov_box <- bxp + 
  stat_pvalue_manual(stat.test, label = "p.adj.signif", tip.length = 0.002, y.position = vsteps, hide.ns = FALSE) +
  stat_n_text(y.pos=-0.001,size = 2.5)
pvcs.sov_box

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p2_legend <- get_legend(pvcs.sov_box)

p3 = grid.arrange(pWT_box + theme(legend.position="none"),
                  psov_box + theme(legend.position="none"),
                  pvcs_box + theme(legend.position="none"),
                  pvcs.sov_box+ theme(legend.position="none"),
                  p2_legend,
                  layout_matrix=rbind(c(1,2,3,4,5)),
                  top=textGrob("mRNA decay rates for uORF genes NMD-targets in wt and decay mutants",gp = gpar(fontsize = 14, fontface = "bold")))

ggsave("~/Desktop/NMD_revise/Figures/Figure 6X. boxplot mRNA decay rates for uORF genes NMD-targets in wt and decay mutants NMDt_rm_TuORFs.pdf",
       plot = p3,
       units = "in",
       width = 13,
       height = 3.5)

