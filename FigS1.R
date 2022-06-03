# Use results from Kallisto quant for CDS regions
library(ggplot2)
library(corrplot)
library(dplyr)

RNAD1 <- read.delim("~/Desktop/kallisto_RNA_quant/D1/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RNAD2 <- read.delim("~/Desktop/kallisto_RNA_quant/D2/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RNAD3 <- read.delim("~/Desktop/kallisto_RNA_quant/D3/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RNAD5 <- read.delim("~/Desktop/kallisto_RNA_quant/D5/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RNAD6 <- read.delim("~/Desktop/kallisto_RNA_quant/D6/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RNAD7 <- read.delim("~/Desktop/kallisto_RNA_quant/D7/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")

RiboD1 <- read.delim("~/Desktop/kallisto_Ribo_quant/D1/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RiboD2 <- read.delim("~/Desktop/kallisto_Ribo_quant/D2/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RiboD3 <- read.delim("~/Desktop/kallisto_Ribo_quant/D3/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RiboD5 <- read.delim("~/Desktop/kallisto_Ribo_quant/D5/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RiboD6 <- read.delim("~/Desktop/kallisto_Ribo_quant/D6/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")
RiboD7 <- read.delim("~/Desktop/kallisto_Ribo_quant/D7/abundance.tsv",header=T,sep="\t",stringsAsFactors = F,quote = "")

TPM=data.frame(RNAD1=RNAD1$tpm,RNAD2=RNAD2$tpm,RNAD3=RNAD3$tpm,RNAD5=RNAD5$tpm,RNAD6=RNAD6$tpm,RNAD7=RNAD7$tpm,
               RiboD1=RiboD1$tpm,RiboD2=RiboD2$tpm,RiboD3=RiboD3$tpm,RiboD5=RiboD5$tpm,RiboD6=RiboD6$tpm,RiboD7=RiboD7$tpm)

TPM$RNA_mean <-rowMeans(TPM[,1:6])
TPM$Ribo_mean <-rowMeans(TPM[,7:12])

TPM2 <- TPM %>% filter(RNA_mean>0.1,Ribo_mean>0.1) %>% filter(RNA_mean<100,Ribo_mean<100)

R2 <- round(cor(TPM2[,1:12]),2)

corrplot(R2, method = "number",type="upper")
