# process Salmon quant data -----------------------------------------------
source("https://bioconductor.org/biocLite.R")
biocLite(c("tximport","GenomicFeatures"))
library(GenomicFeatures)
library(tximport)
library(readr)
library(DESeq2)
library(biomaRt)

samples = read.csv("./Data/RAW/ASD_brain/samples.csv", row.names = 1)
files = file.path("./Data/RAW/ASD_brain", paste0(rownames(samples), "_quant"), "quant.sf")
names(files) = samples$Names

txdb = makeTxDbFromGFF("./Data/RAW/Annotation/gencode.v25.annotation.gtf.gz")
k = keys(txdb, keytype = "GENEID")
df = select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene = df[ ,2:1]
head(tx2gene)

txi.salmon = tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)
colnames(txi.salmon$counts)= rownames(samples)

dds = DESeqDataSetFromTximport(txi.salmon, samples, ~Class)
dds_DE =DESeq(dds)
res05 =results(dds_DE, alpha = 0.05, contrast = c("Class","ASD","Con")) 

pdf("./Figures/MA_plot.pdf")
plotMA(res05, main="MA plot of ASD fold change", ylim=c(-3,3))
dev.off()
res05= res05[!is.na(res05$padj), ]  
rownames(res05) = unlist(lapply(strsplit(rownames(res05), "\\."),"[[",1)) ## remove trailing decimals

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","entrezgene","gene_biotype","start_position","end_position","chromosome_name", "external_gene_name"),
            filters = "ensembl_gene_id", values=rownames(res05))

sig_match = match(rownames(res05),  map$ensembl_gene_id)
res05$gene_name = map$external_gene_name[sig_match]
res05$biotype = map$gene_biotype[sig_match]
res05$chromosome = map$chromosome_name[sig_match]
res05$start_pos = map$start_position[sig_match]
res05$end_pos = map$end_position[sig_match]
res05$entrez = map$entrezgene[sig_match]

lncRNA_filter = c("3prime_overlapping_ncrna", "antisense","antisense RNA", "lincRNA","ncrna host","processed_transcript", "sense_intronic" , "sense_overlapping")

lnc_test=c()
for (biotype in res05$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0 ) )
}
res05$lncRNA = lnc_test
table(res05$lncRNA)



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

res05 =res05[res05$padj <= 0.05, ]
