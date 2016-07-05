setwd("C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/DESeq2")
list.files()

library(dplyr)

DEG = read.csv("DEG_Results_alpha05.csv")


lncRNAs = DEG$lncRNA == TRUE

table( DEG$biotype[lncRNAs])

lnc_types = table(DEG$biotype[lncRNAs])
lnc_types = lnc_types[lnc_types != 0]

## 293 lncRNAs / 1692 PC-genes

pdf("LncRNA Biotypes.pdf")
par(mar = c(12, 4, 4, 2));
barplot(lnc_types, las=2, main = "Biotype of (293) LncRNAs Differentially Expressed in ASD")
dev.off()


lnc_xp = group_by(DEG[lncRNAs,], DEG$biotype[lncRNAs])
lnc_types = lnc_types[lnc_types != 0]

summarise(lnc_xp, mean(log2FoldChange))

summarise(lnc_xp, mean(log2FoldChange))
# Source: local data frame [6 x 2]

# DEG$biotype[lncRNAs] mean(log2FoldChange)
# (fctr)                (dbl)
# 1 3prime_overlapping_ncrna            1.7473859
# 2                antisense           -0.3338986
# 3                  lincRNA           -0.2710390
# 4     processed_transcript           -0.2203220
# 5           sense_intronic           -0.5102447
# 6        sense_overlapping            0.5342363

# quantile(lnc_xp$log2FoldChange)
# 0%       25%       50%       75%      100% 
# -3.308318 -1.689759 -1.305628  1.533424  3.155564 

# table(lnc_xp$log2FoldChange < 0)
# FALSE  TRUE 
# 124   169 

# 169/ (169+124)
# [1] 0.5767918 % down-regulated


#### Protein-coding genes
quantile(DEG$log2FoldChange[!lncRNAs])
# 0%        25%        50%        75%       100% 
# -4.2818118 -1.4451784  0.8471417  1.3827251  4.4462584 

source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Hsapiens.UCSC.hg38.knownGene")

library(TxDb.Hsapiens.UCSC.hg38.knownGene)
genes = genes(TxDb.Hsapiens.UCSC.hg38.knownGene)

lncs = DEG[lncRNAs, c(10,11,12)]
lncs$chromosome = paste0("chr", lncs$chromosome)

gr = makeGRangesFromDataFrame(lncs, start.field = "start_pos", end.field = "end_pos")


nearby = genes[nearest(gr, genes)]

library(biomaRt)
mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","gene_biotype","chromosome_name","start_position","end_position"),
            filters = "entrezgene", values =  nearby$gene_id)
### makes genes from Bspan list?
