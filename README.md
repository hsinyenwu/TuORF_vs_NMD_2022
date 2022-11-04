### The code provided here is for the following paper:  
[Nonsense-mediated decay is not a major mechanism for regulating the uORF-containing mRNAs in Arabidopsis](https://www.biorxiv.org/content/10.1101/2021.09.16.460672v3)

**We provided six lines of evidence that NMD is not a major mechanism regulating TuORF mRNAs in plants:**  
(1) NMD targets, but not TuORF mRNAs, are upregulated in various NMD mutants (Figure 3).  
(2) TuORF mRNAs have significantly higher transcript levels, while the NMD targets have lower transcript levels than other genes (Figure 4A and 4E).  
(3) NMD and TuORFs additively repress translation, suggesting they repress translation through separate mechanisms (Figure 4B, 4F).  
(4) While NMD targets are efficiently degraded through both 5’ to 3’ and 3’ to 5’ decay, TuORF mRNAs, like other genes, mainly rely on 5’ to 3’ decay (Figure 5).  
(5) The decay rates of NMD targets are consistently higher than those of TuORF mRNAs (Figure S11).  
(6) NMD targets also appear to be degraded by additional pathway(s) in the absence of both VCS (5'->3') and SOV (3'->5') pathways (Figure S11H).  

The data has been uploaded to [GEO](https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&cad=rja&uact=8&ved=2ahUKEwiUjJLUz5n6AhVrjokEHdoyCFIQFnoECAQQAQ&url=https%3A%2F%2Fwww.ncbi.nlm.nih.gov%2Fgeo%2F&usg=AOvVaw3Dc0qZ8-sNL7LwyPCWsoAr) under GSE183264. The datasets will be released to public after the paper is accepted. 
The Arabidopsis annotation Araport11 and TAIR10 genome sequence files can be found on [TAIR](https://www.arabidopsis.org)

The TuORF positions can be visualized on the [JBrowse TAIR](https://jbrowse.arabidopsis.org/index.html?data=Araport11&loc=Chr1%3A21537..32850&tracks=TAIR10_genome%2CA11-GL-Jul22%2CA11-PC-Jul22%2CSALK_tDNAs&highlight=) after the paper is accepted.

The downstream analysis is done in R and used numorous R packages. Please see each file for detail. For visualizing the distribution of Ribo-seq reads for individual genes, please use [RiboPlotR](https://github.com/hsinyenwu/RiboPlotR).  


#### FIGURES  
Figure 1. Enhanced Ribo-seq coverage improves uORF identification  
Figure 2. TuORF-containing genes have higher mRNA abundance but are associated with lower translation efficiency, shorter mRNA half-lives, and lower protein abundance  
Figure 3: TuORF mRNAs are not upregulated in the NMD mutants  
Figure 4. Expression patterns of NMD targets and TuORF mRNAs  
Figure 5. TuORF mRNAs and NMD targets are degraded through different mRNA decay mechanisms  

#### SUPPLEMENTARY DATA  
Figure S1. Correlations among samples and comparisons between our previous and current datasets  
Figure S2. Examples of TuORFs in important regulatory genes  
Figure S3. Protein abundance of the uORF-containing genes in Arabidopsis root  
Figure S4. The mRNA levels and translation efficiency of TuORF mRNAs in tomato  
Figure S5. mRNA levels of uORF-containing genes in different ecotypes, growth stages, and tissues in Arabidopsis and tomato  
Figure S6. mRNA levels of CPuORF-containing genes compared to those of no-uORF genes  
Figure S7. (A) TuORF translation levels are positively correlated with the mRNA levels. (B) The distributions of the mRNA levels of TuORF and UuORF genes in Arabidopsis seedlings. Vertical dashed lines indicate the median value of each group  
Figure S8. TuORF genes are larger than other genes  
Figure S9. Lengths of the uORF peptides for NMD target and non-NMD target TuORFs  
Figure S10. Lengths of UTRs and CDS for NMD targets and TuORF genes  
Figure S11. The mRNA decay rates of NMD targets and TuORF-containing mRNAs in different genetic backgrounds  

#### TABLES
Table S1. Translated ORFs identified by RiboTaper  
•	Table S1A-C. Translated uORFs identified in the current dataset and our previous Arabidopsis shoot and root datasets  
•	Table S1D-F. Translated mORFs identified in the current dataset and our previous Arabidopsis shoot and root datasets  
Table S2. GO-term analysis of TuORF genes  
Table S3. Summary of the large-scale datasets used in this study
Table S4. Expression data TPM values for RNA-seq, Ribo-seq and TE for each gene (CDS region)
Table S5. NMD target gene list  
•	Table S5A. High-confidence NMD targets expressed in our data  
•	Table S5B. 49 genes identified as high-confidence NMD target TuORF genes  

For more about Ribo-seq library construction, please see [here](https://github.com/hsinyenwu/Riboseq_protocol_2022).  
For more about Ribo-seq visualization, please see [here](https://github.com/hsinyenwu/RiboPlotR).  

#### Files used in this analysis:
Table S4: TE_CDS data, will be released after paper accepted for publication  
Annotation file used here: [Araport11_20181206_max_isoform.gtf](https://github.com/hsinyenwu/TuORF_vs_NMD_2022/blob/main/Data/Araport11_20181206_max_isoform.gtf.zip)  
Arabidopsis expression data:
[E-MTAB-7978](https://www.ebi.ac.uk/gxa/experiments/E-MTAB-7978/Downloads)




