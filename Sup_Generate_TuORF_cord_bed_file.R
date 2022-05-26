```
rm(list=ls())
library(GenomicRanges)
library(GenomicFeatures)
library(Rsamtools)
library(ORFik)
library(rtracklayer)
library(dplyr)

# Load genome sequences
FA <- FaFile("~/Desktop/Leaky_scanning/TAIR10_chr_all_2.fas")
# Generate txdb object
txdb <- makeTxDbFromGFF("~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
# Load GTF file
GTF <- read.delim("~/Desktop/CTRL_v1/Araport11+CTRL_20181206.gtf",header=F,stringsAsFactors = F)
# Generate 5'UTR ranges
fiveUTR <- fiveUTRsByTranscript(txdb, use.names=T)
# Generate CDS ranges
CDS <- cdsBy(txdb,by="tx", use.names=T)

#load RiboTaper output ORFs_max_filt file

Exp <- read.delim(file="~/Desktop/uORFs_miRNA/ORFs_max_filt_Araport_both_round_reads",header=T,stringsAsFactors=F,sep="\t")
Exp_uORFs <- Exp %>% filter(category=="uORF") %>% mutate(ORF_pept_length=nchar(ORF_pept))
nrow(Exp_uORFs)

addCdsPhase <- function(cds_by_tx)
{
  cds_phase <- pc(rep(IntegerList(0), length(cds_by_tx)),
                  heads(cumsum(width(cds_by_tx)) %% 3L, n=-1L))
  unlisted_cds_by_tx <- unlist(cds_by_tx, use.names=FALSE)
  mcols(unlisted_cds_by_tx)$phase <- unlist(cds_phase, use.names=FALSE)
  relist(unlisted_cds_by_tx, cds_by_tx)
}

for(i in 1:nrow(Exp_uORFs)){
  #if(i%%100==0) print(i)
  print(i)
  gc()
  # a is the 5'UTR range for the ith uORF
  a=fiveUTR[Exp_uORFs$transcript_id[i]]
  # b is the CDS range for the ith uORF
  b=CDS[Exp_uORFs$transcript_id[i]]
  # d is the output ranges for the findUORFs function
  d=findUORFs(fiveUTRs=a,fa=FA, startCodon="ATG",longestORF=F,cds=b,restrictUpstreamToTx=T)
  # then identify the RiboTaper identified uORFs that is identical with which findUORFs defined uORFs 
  Pept <- as.character(translate(extractTranscriptSeqs(FA,d)))
  Pept_remove_stop <- gsub("\\*","",Pept)
  Pept_uORF <- d[which(Pept_remove_stop==Exp_uORFs[i,]$ORF_pept)]
  # Make the mRNA rows for gtf file
  mRNAr = GTF %>% filter(V3=="mRNA",grepl(pattern=Exp_uORFs$transcript_id[i],V9))
  # Make the CDS rows for gtf file
  CDSr <- export.gff(Pept_uORF,con=paste0("~/Desktop/uORF_gtfs/cds.gtf"),format="gff2")
  CDSr <- read.delim("~/Desktop/uORF_gtfs/cds.gtf",skip=4,header=F)
  CDSr$V3 <- "CDS"
  CDSr$V9 <- paste0("gene_id ","\"",substr(Exp_uORFs$transcript_id[i],1,9),"\"; ", "transcript_id ","\"",Exp_uORFs$transcript_id[i],"\";")
  #output the file
  write.table(rbind(mRNAr,CDSr),paste0("~/Desktop/uORF_gtfs/","uORF_temp",".gtf"),quote = F, col.names = F, row.names = F, sep = "\t",append = F)
  txdb2 <- makeTxDbFromGFF("~/Desktop/uORF_gtfs/uORF_temp.gtf",format="gtf", dataSource="Araport11",organism="Arabidopsis")
  CDS2 <- cdsBy(txdb2,by="tx", use.names=T)
  CDS2 <- addCdsPhase(CDS2)
  export.bed(CDS2,paste0("~/Desktop/uORF_gtfs/BED/",Exp_uORFs$ORF_id_tr[i],".bed"))
  BED <- read.delim(file=paste0("~/Desktop/uORF_gtfs/BED/",Exp_uORFs$ORF_id_tr[i],".bed"),header=F,stringsAsFactors=F,sep="\t")
  BED$V4 <- Exp_uORFs$ORF_id_tr[i]
  write.table(BED,"~/Desktop/uORF_gtfs/TuORFs.bed",quote = F, col.names = F, row.names = F, sep = "\t",append = T)
}

```
