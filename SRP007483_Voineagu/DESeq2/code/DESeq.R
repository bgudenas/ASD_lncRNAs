setwd("C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/DESeq2")

library(DESeq2)
library(biomaRt)


counts=read.csv("C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/RAW/Counts_table.csv",row.names=1)
colnames(counts)=c("ASD_3334","ASD_3536","ASD_3738","CON_3940","CON_4142","CON_4344")

counts=counts[1:(nrow(counts)-5),]   ## remove bottom 5 quality assurance rows from HTseq

condition=data.frame(condition=c("ASD","ASD","ASD","CON","CON","CON"))
rownames(condition)=colnames(counts)

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","entrezgene","gene_biotype","start_position","end_position","chromosome_name"),
            filters = "ensembl_gene_id", values=rownames(counts))

dds=DESeqDataSetFromMatrix(countData = counts ,colData = condition, design = ~ condition)
dds_DE =DESeq(dds)
res05 =results(dds_DE, alpha = 0.05)
# sum(res05$padj < 0.05, na.rm=TRUE)
# [1] 1985

plotMA(res05, main="DESeq2")


res05 = res05[!is.na(res05$padj),]   ### remove genes with NA as p-adjusted
res_sig = res05[ res05$padj <= 0.05,]  ## filter out genes above alpha 0.05

sig_match = match(rownames(res_sig),  map$ensembl_gene_id)

map$gene_biotype[sig_match]

res_sig$biotype = map$gene_biotype[sig_match]

### lncRNA_filter represents all possible lncRNA biotypes as defined by Ensembl
lncRNA_filter = c("3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript", "sense_intronic" , "sense_overlapping")

lnc_test=c()
for (biotype in res_sig$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0 ) )
}
res_sig$lncRNA = lnc_test

res_sig$chromosome = map$chromosome_name[sig_match]
res_sig$start_pos = map$start_position[sig_match]
res_sig$end_pos = map$end_position[sig_match]
res_sig$entrez = map$entrezgene[sig_match]

resOrdered <- res_sig[order(res_sig$padj),]

##----- write out significant DEG
write.csv(resOrdered,"DEG_Results_alpha05.csv")


rld = rlog(dds_DE)
sampleDists = dist(t(assay(rld)))

h_rld = hclust(sampleDists)
plot(h_rld, main="Clustering of samples based on Rlog gene-count")

library("RColorBrewer")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("pheatmap")
library("pheatmap")

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(dds)
colnames(sampleDistMatrix) <- NULL

colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
