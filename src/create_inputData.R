
# Create Supp data file ---------------------------------------------------

library(GenomicFeatures)
library(tximport)
library(readr)
DEG_samples = read.csv("./Data/RAW/ASD_brain/samples.csv", row.names = 1) ##Sample metadata w/ rownames as sample names

## filepath to quantification files produced by Salmon
files = file.path("./Data/RAW/ASD_brain", paste0(rownames(DEG_samples), "_quant"), "quant.sf")
names(files) = rownames(samples)

## Gene Annotation file from  ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/gencode.v25.annotation.gtf.gz
txdb = makeTxDbFromGFF("./Data/RAW/Annotation/gencode.v25.annotation.gtf.gz") 

k = keys(txdb, keytype = "GENEID")
df = select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene = df[ ,2:1]
head(tx2gene) #transcript to gene map

txi.salmon = tximport(files, type = "salmon", tx2gene = tx2gene, reader = read_tsv)

# Add in Brainspan Data ---------------------------------------------------

Bspan_rows<-read.csv("./Data/RAW/Brainspan/rows_metadata.csv")  ##genelist of Brainspan transcriptome summarized to genes,5/15
SFARI <-read.csv("./Data/RAW/SFARI/gene-summary.csv") #7/11/15 SFARI ASD genes, 707 total
SFARI_scores<-read.csv("./Data/RAW/SFARI/gene-score.csv") #SFARI ASD risk gene scores


## load expression matrix and sample metadata
Expr <- read.csv("./Data/Raw/Brainspan/expression_matrix.csv",header=FALSE,row.names=1)  ##RPKM expression data
clinical <- read.csv("./Data/Raw/Brainspan/columns_metadata.csv",row.names=1)  ##matrix containing sample metadata
ME16 = read.csv("./Data/RAW/parikshak_m16.csv")
ME16 = ME16$ENSEMBL.GENE.ID

GTEx = read.csv("./Data/RAW/GTEx/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.csv")
# GTEx_samples from http://www.gtexportal.org/home/tissueSummaryPage
GTEx_samples = read.csv("./Data/RAW/GTEx//GTEx_sample_metadata.csv")

CNVs = read.csv("./Data/RAW/SFARI/cnv-summary.csv", row.names = NULL)
colnames(CNVs)=c("Chromosome","CNV.Locus","CNV.Type","Largest.CNV.Size","PMID" ,"Title","Author","Year" ,"Report.Class")
CNVs = CNVs[CNVs$Report.Class == "Major", ]

rm(txdb, k,files,first_time,df,tx2gene)
save.image("./Data/All_Sup_Data.RData")
