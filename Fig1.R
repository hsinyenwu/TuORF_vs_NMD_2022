#############################################################################################
# Figure 1A meta plot for 3 nt periodicity
#############################################################################################
# Use the metaplot function from RiboTaper

#############################################################################################
# Figure 1B and 1C
#############################################################################################
# Performed using RiboseqQC (https://github.com/ohlerlab/RiboseQC)

#############################################################################################
# Figure 1D read diversity
#############################################################################################

# To determine the diversity of read counts in our new data versus the root and shoot data from Hsu et al. 2016, 
# we randomly select 1 to 10 million ribo-seq reads mapped to chromosome 1 in the three datasets and identified the number of unique P-sites. 
# The more unique P-sites suggest more diverse ribosome footprints.

rm(list=ls())
library(dplyr)
Root <- read.delim("~/Desktop/BACKUP/Root_P_sites_all_4Col.bed",header=F,sep="\t",stringsAsFactors = F,skip=0) #2,157,859 rows, sum total reads 132,308,885
Shoot <- read.delim("~/Desktop/BACKUP/Shoot_P_sites_all_4Col.bed",header=F,sep="\t",stringsAsFactors = F,skip=0) #1,908,790 rows, sum total reads 106,752,951
NewCTRL <- read.delim("/Users/wu/Desktop/CTRL_v1/CTRL_TPM0.25_P_sites_sort_count",header=F,sep="\t",stringsAsFactors = F,skip=0) #11,209,299 rows, sum total reads 298,004,424

RootChr1 <- Root %>% filter(V2==1)
nrow(RootChr1) #549847
ShootChr1 <- Shoot %>% filter(V2==1)
nrow(ShootChr1) #485377
NewCTRLChr1 <- NewCTRL %>% filter(V2==1)
nrow(NewCTRLChr1) #2913663

RootChr1R <- as.data.frame(lapply(RootChr1, rep, RootChr1$V1))
nrow(RootChr1R) #31048876
ShootChr1R <- as.data.frame(lapply(ShootChr1, rep, ShootChr1$V1))
nrow(ShootChr1R) #26708668
NewCTRLChr1R <- as.data.frame(lapply(NewCTRLChr1, rep, NewCTRLChr1$V1))
nrow(NewCTRLChr1R) #73870022

NUM <- seq(1e+6,2e+7, by = 1e+6)

RootChr1Cdf <- c()
for(i in 1:20) {
  RootChr1C <- RootChr1R[sample(nrow(RootChr1R), NUM[i]), ]
  res <- RootChr1C[!duplicated(RootChr1C), ]
  RootChr1Cdf[i] <- nrow(res)
}
RootChr1Cdf

ShootChr1Cdf <- c()
for(i in 1:20) {
  ShootChr1C <- ShootChr1R[sample(nrow(ShootChr1R), NUM[i]), ]
  res <- ShootChr1C[!duplicated(ShootChr1C), ]
  ShootChr1Cdf[i] <- nrow(res)
}
ShootChr1Cdf

NewCTRLChr1Cdf <- c()
for(i in 1:20) {
  NewCTRLChr1C <- NewCTRLChr1R[sample(nrow(NewCTRLChr1R), NUM[i]), ]
  res <- NewCTRLChr1C[!duplicated(NewCTRLChr1C), ]
  NewCTRLChr1Cdf[i] <- nrow(res)
}
NewCTRLChr1Cdf

dfCTRL <- data.frame(total_read=NUM, distinct_p_site=NewCTRLChr1Cdf,Samples="At_CTRL")
dfRoot <- data.frame(total_read=NUM, distinct_p_site=RootChr1Cdf,Samples="At_Root")
dfShoot <- data.frame(total_read=NUM, distinct_p_site=ShootChr1Cdf,Samples="At_Shoot")

df <- rbind(dfCTRL,dfRoot,dfShoot)
save(df,"~/Desktop/atRTD3/readDiversity.RData")
df$distinct_p_site <- df$distinct_p_site/1000000
p3 <- ggplot(df, aes(x=total_read,y=distinct_p_site, colour = Samples)) +
  geom_line(aes(color=Samples),size=0.7) +
  geom_point(aes(color=Samples)) +
  theme_bw() +
  labs(x="Million footprint counts",y="Million distinct P-sites in 20M samples")
p3
ggsave("~/Desktop/atRTD3/Footprint_diversity3.pdf") 

#############################################################################################
# Figure 1E ORF identified
#############################################################################################


#############################################################################################
# Figure 1F,1G uORF examples
#############################################################################################

library(RiboPlotR)
gene.structure(annotation="~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")

rna_bam.ribo(RNAseqBam1="~/Desktop/CTRL_v1/RNA_CTRL_merged.bam",
             Ribo1="/Users/wu/Desktop/CTRL_v1/CTRL_expressed_P_sites_sort_count",
             RNAlab1="Ctrl_RNA count",
             Ribolab1="Ctrl_ribo count",
             S_NAME1="Ctrl",
             RNAbackground="#FFFFE0",
             RNAseqBamPaired="paired")
PLOTc("AT1G05160",isoform=2,NAME = "ATKAO1")
PLOTc("AT3G02470",isoform=3,NAME = "SAMDC1")

# use the following to plot uORF regions
# uorf.structure(uorf_annotation= "~/Desktop/uORF.gtf", format="gtf",dataSource="Araport",organism="Arabidopsis thaliana")
# PLOTc("AT1G05160", isoform=2, uORF="AT1G05160", NAME = "ATKAO1")
# PLOTc("AT3G02470", isoform=3, uORF="AT3G02470", NAME = "SAMDC1")

#gtf for AT1G05160 uORF
#1	Araport11	mRNA	1487378	1491005	.	-	.	gene_id AT1G05160; transcript_id AT1G05160.1; gene_biotype protein_coding;
#1	rtracklayer	CDS	1490730	1490837	.	-	.	gene_id "AT1G05160"; transcript_id "AT1G05160.1";
#1	rtracklayer	CDS	1489912	1489953	.	-	.	gene_id "AT1G05160"; transcript_id "AT1G05160.1";

#gtf for AT3G02470 uORF
#3	Araport11	mRNA	509096	511489	.	+	.	gene_id AT3G02470; transcript_id AT3G02470.3; gene_biotype protein_coding;
#3	rtracklayer	CDS	509815	509952	.	+	.	gene_id "AT3G02470"; transcript_id "AT3G02470.3";
#3	rtracklayer	CDS	510048	510068	.	+	.	gene_id "AT3G02470"; transcript_id "AT3G02470.3";





