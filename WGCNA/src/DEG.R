
# Take HTSeq output and concatenate a counts table ------------------------

#path to main project directory
root = "C:/Users/Brian/Google Drive/ASD_lncRNAs/WGCNA"

### change workdir to where raw count data is
setwd(paste0(root, "/RAW/SRP007483/"))

count=c()
samples=c()
for (file in list.files()){
    sample=strsplit(strsplit(file,"_")[[1]][3],".txt")[[1]]
    status = strsplit(file, "_")[[1]][2]
    sample = paste(status, sample, sep="_")
    
    samples=c(samples,sample)
    
    counts=read.table(file)
    genes=counts$V1
    count=c(count,counts$V2)
}

counts = matrix(data=count,nrow=length(genes),ncol=length(list.files()))
rownames(counts)=genes
colnames(counts)=samples

setwd(root)

# Differential Expression using DESeq2 ------------------------------------
library(DESeq2)
library(biomaRt)



counts = counts[1:(nrow(counts)-5), ]   ## remove bottom 5 quality assurance rows produced by HTseq

condition=data.frame(condition=c("ASD","ASD","ASD","CON","CON","CON"))
rownames(condition)=colnames(counts)


mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","entrezgene","gene_biotype","start_position","end_position","chromosome_name", "external_gene_name"),
            filters = "ensembl_gene_id", values=rownames(counts))

dds=DESeqDataSetFromMatrix(countData = counts ,colData = condition, design = ~ condition)
dds_DE =DESeq(dds)
res05 =results(dds_DE, alpha = 0.05, contrast = c("condition","ASD","CON")) ## means fold change = ASD/CON

# MA plot to test normalization  ------------------------------------------
pdf("./Figures/MA_plot.pdf")
plotMA(res05, main="MA plot of ASD fold change", ylim=c(-3,3))
dev.off()

res05= res05[!is.na(res05$padj), ]   ### remove genes with NA as p-adjusted

sig_match = match(rownames(res05),  map$ensembl_gene_id)

res05$gene_name = map$external_gene_name[sig_match]
res05$biotype = map$gene_biotype[sig_match]
res05$chromosome = map$chromosome_name[sig_match]
res05$start_pos = map$start_position[sig_match]
res05$end_pos = map$end_position[sig_match]
res05$entrez = map$entrezgene[sig_match]
### lncRNA_filter represents all possible lncRNA biotypes as defined by Ensembl
lncRNA_filter = c("3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript", "sense_intronic" , "sense_overlapping")

lnc_test=c()
for (biotype in res05$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0 ) )
}
res05$lncRNA = lnc_test


# Make Volcano plot for DEG visualization ---------------------------------
tab = data.frame(logFC = res05$log2FoldChange, negLogPval = -log10(res05$padj))

pdf("./Figures/Volcano.pdf")

par(mar = c(5, 4, 4, 4))
plot(tab, cex = 0.6, xlab=expression(Log[2]~fold~change), pch=1,
                ylab = expression(-Log[10]~pvalue), ylim = c(0,13))
title(main = "Differentially Expressed Genes in the ASD brain")
lfc = 1
pval = 0.05
signGenes = (abs(tab$logFC) >= lfc & tab$negLogPval > -log10(pval) & !res05$lncRNA)
points(tab[signGenes, ], pch = 1, cex = 0.8, col = "red3")
signLncRNAs = (abs(tab$logFC) >= lfc & tab$negLogPval > -log10(pval) & res05$lncRNA)
points(tab[signLncRNAs, ], pch = 1, cex = 0.8, col = "blue")
abline(h = -log10(pval), col = "magenta", lty = 2, lwd = 2)
abline(v = c(-lfc, lfc), col = "green2", lty = 2, lwd= 2)
mtext(paste("pval =", pval), side = 4, at = -log10(pval), cex = 0.8, line = 0.5, las = 1)
mtext(c(paste("-", lfc, "fold"), paste("+", lfc, "fold")), side = 3, at = c(-lfc, lfc),
      cex = 0.8, line = 0.5)

dev.off()



# Filter out non-significant genes ----------------------------------------
res_sig = res05[ res05$padj <= 0.05 & abs(res05$log2FoldChange) >=1  ,] 


res_sig <- res_sig[order(res_sig$padj),]

### check for genes not found in newest version of ensembl
#   table(is.na(res_sig$biotype))
#   FALSE  TRUE 
#   1980     5

# remove genes not found in newest ensembl --------------------------------
res_sig = res_sig[!is.na(res_sig$biotype),  ]

#     nrow(res_sig)
#     [1] 1862
#     table(res_sig$lncRNA)
#     FALSE  TRUE 
#     1584   278



##----- write out significant DEG
write.csv(res_sig,"./Data/DEG_Results_alpha05.csv")


# check for sample outliers by clustering ---------------------------------

rld = rlog(dds_DE)
sampleDists = dist(t(assay(rld)))

h_rld = hclust(sampleDists, method="average")
plot(h_rld, main="Clustering of samples based on Rlog gene-count")
## no discernible outliers

