#Figure 5 dot plot with error bars
library(openxlsx)
library(ggplot2)

############
# April 23

DLA1 <- read.xlsx("~/Desktop/NMD_revise/DLA_results_Phong_05062022.xlsx",sheet = 1)

DLA1$group <- as.factor(paste(DLA1$Line, DLA1$Plasmid))
DLA1$group <- factor(DLA1$group, levels = c("Col AUG", "Col AGG", "lba1 AUG","lba1 AGG"))

p1 <- ggplot(DLA1, aes(x=group, y=Ratio,color=group)) + 
  geom_boxplot() +theme_classic() + ggtitle("AT2G01570 (RGA1) Apr23")+labs(tag = "A") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")+ylim(c(0,9))

p1

DLA2 <- read.xlsx("~/Desktop/NMD_revise/DLA_results_Phong_05062022.xlsx",sheet = 2)

DLA2$group <- paste(DLA2$Line, DLA2$Plasmid)
DLA2$group <- factor(DLA2$group, levels = c("Col AUG", "Col AGG", "lba1 AUG","lba1 AGG"))

p2 <- ggplot(DLA2, aes(x=group, y=Ratio,color=group)) + 
  geom_boxplot() +theme_classic() + ggtitle("AT3G08730 (S6K1) Apr23") +labs(tag = "C")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")+ylim(c(0,7))

p2
############
# May 6

DLA3 <- read.xlsx("~/Desktop/NMD_revise/DLA_results_Phong_05062022.xlsx",sheet = 3)

DLA3$group <- as.factor(paste(DLA3$Line, DLA3$Plasmid))
DLA3$group <- factor(DLA3$group, levels = c("Col AUG", "Col AGG", "lba1 AUG","lba1 AGG"))

p3 <- ggplot(DLA3, aes(x=group, y=Ratio,color=group)) + 
  geom_boxplot() +theme_classic() + ggtitle("AT2G01570 (RGA1) May6")+labs(tag = "B") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")+ylim(c(0,2))

p3



DLA4 <- read.xlsx("~/Desktop/NMD_revise/DLA_results_Phong_05062022.xlsx",sheet = 4)

DLA4$group <- paste(DLA4$Line, DLA4$Plasmid)
DLA4$group <- factor(DLA4$group, levels = c("Col AUG", "Col AGG", "lba1 AUG","lba1 AGG"))

p4 <- ggplot(DLA4, aes(x=group, y=Ratio,color=group)) + 
  geom_boxplot() +theme_classic() + ggtitle("AT3G08730 (S6K1) May6") +labs(tag = "D") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + xlab("")+ylim(c(0,1.5))

p4

get_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

p4_legend <- get_legend(p4)

p5 <- grid.arrange(p1 + theme(legend.position="none"),
                   p3 + theme(legend.position="none"),
                   p2 + theme(legend.position="none"),
                   p4 + theme(legend.position="none"),
                   p4_legend,
                   layout_matrix=rbind(c(1,2,3,4,5)),
                   top=textGrob("   Figure 5 ", gp=gpar(fontsize=13), x = 0, hjust = 0))

ggsave("~/Desktop/NMD_revise/Figure 5x DLA.pdf",
       plot = p5,
       units = "in",
       width = 14.5,
       height = 5)
