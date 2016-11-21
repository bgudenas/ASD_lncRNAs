#9/20/16 written by Brian Gudenas
# WGCNA analysis

library(WGCNA)
library(stringr)
library(biomaRt)
library(GenomicRanges)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

setwd("./WGCNA")

Bspan_rows<-read.csv("./RAW/Brainspan/rows_metadata.csv")  ##genelist of Brainspan transcriptome summarized to genes,5/15
DEG <-read.csv("./Data/DEG_Results_alpha05.csv",row.names = 1) ## ASD Differentially expression 
SFARI <-read.csv("./RAW/SFARI/gene-summary.csv") #7/11/15 SFARI ASD genes, 707 total
scores<-read.csv("./RAW/SFARI/gene-score.csv") #SFARI ASD risk gene scores

score_match = match(scores$Gene.Symbol, SFARI$Gene.Symbol)
scores$entrez = SFARI$Entrez.GeneID[score_match]
table(DEG$biotype[DEG$lncRNA])/sum(DEG$lncRNA)*100

## load expression matrix and sample metadata
Expr <- read.csv("./Raw/Brainspan/expression_matrix.csv",header=FALSE,row.names=1)  ##RPKM expression data
clinical <- read.csv("./Raw/Brainspan/columns_metadata.csv",row.names=1)  ##matrix containing sample metadata



# convert age based on diff units to a single continous variable -- Months Post-conception (conception ~10 months)
clinical$Months.Post.Conception = 0
for (i in 1:nrow(clinical)) {
    age = as.numeric(strsplit(clinical$age[i], " ")[[1]][1])
    unit = strsplit(clinical$age[i], " ")[[1]][2]
    
    if (unit == "pcw"){ age = age/4 #convert to months
        }else if (unit == "mos") { age = age+10 ## already in months but need to add 10 months for conception
        }else if (unit == "yrs") { age = (age*12)+10}## convert years to months then add 10 for conception
    clinical$Months.Post.Conception[i]=age
}

# create unique IDs for Expr and clinical samples -------------------------
colnames(Expr)<-str_c(clinical$structure_acronym, clinical$Months.Post.Conception, clinical$gender, clinical$donor_id, sep="_")
rownames(clinical) = colnames(Expr)

# filter samples to only retain regions within the neocortex by matching " cortex" and not cerebellar cortex which is in the cerebellum and not neocortex
cortical_regions = grepl(" cortex", clinical$structure_name) & clinical$structure_name != "cerebellar cortex"
# double-check brain regions included in analysis
table(clinical$structure_name[cortical_regions])
#Remove the  primary motor-sensory cortex bc has only 5 samples while the rest all have ~30
cortical_regions = grepl(" cortex", clinical$structure_name) & clinical$structure_name != "cerebellar cortex" & clinical$structure_name != "primary motor-sensory cortex (samples)"

## total number of samples included for downstream analysis
sum(cortical_regions)
# [1] 352

# filter clinical and Expr by cortical_region -----------------------------
clinical = clinical[cortical_regions, ]

Expr = Expr[ ,cortical_regions] 


datExpr0 = t(Expr)
colnames(datExpr0)=Bspan_rows$ensembl_gene_id

#-- THis loop checks the # of times each gene is expressed in a sample >= 1 RPKM
vals=c()
for (i in 1:ncol(datExpr0)) {
    vals=c(vals,sum(datExpr0[,i] >= 1) )
}



Gene_filter = (vals >= 3 )  ## genes required to have normalized RPKM >=1 in atleast 3/352 of samples
##number of kept genes
# [1] 23668

datExpr0 = datExpr0[, Gene_filter]


##--Create Genelist dataframe containing all gene info from (DEG and SFARI)
genelist = Bspan_rows[Gene_filter,]
DEG_map = match(genelist$ensembl_gene_id, rownames(DEG))


genelist$L2FC = DEG$log2FoldChange[DEG_map]
genelist$L2FC[is.na(genelist$L2FC)] = 0    ### Give 0 to any gene not Differentially expressed

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","band","gene_biotype","chromosome_name","start_position","end_position"),
            filters = "ensembl_gene_id", values =  genelist$ensembl_gene_id)


mart_map = match(genelist$ensembl_gene_id, map$ensembl_gene_id)

genelist$biotype = map$gene_biotype[mart_map]
genelist$chromosome = map$chromosome_name[mart_map]
genelist$start = map$start_position[mart_map]
genelist$end = map$end_position[mart_map]

ASD_match = match(genelist$entrez_id, scores$entrez)
genelist$ASD_score = scores$Score[ASD_match]


lncRNA_filter = c("3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript", "sense_intronic" , "sense_overlapping")

lnc_test=c()
for (biotype in genelist$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0))
}
genelist$lncRNA = lnc_test


LncRNA_Gr = makeGRangesFromDataFrame(DEG[DEG$lncRNA, colnames(DEG)=="chromosome" | colnames(DEG)=="start_pos" |  colnames(DEG)=="end_pos"], start.field = "start_pos", end.field = "end_pos" )
genelist_Gr = makeGRangesFromDataFrame(genelist[genelist$lncRNA != TRUE & !is.na(genelist$start), colnames(genelist)=="chromosome" | colnames(genelist)=="start" |  colnames(genelist)=="end" ], start.field = "start", end.field = "end" )
names(genelist_Gr) = genelist$gene_symbol[genelist$lncRNA != TRUE & !is.na(genelist$start) ]

nearest(x = LncRNA_Gr, subject = genelist_Gr)
near_genes = names(genelist_Gr)[nearest(x = LncRNA_Gr, subject = genelist_Gr)]
DEG$nearest[DEG$lncRNA] = near_genes
DEG$nearest_ASD_score[DEG$lncRNA] = genelist$ASD_score[match( near_genes, genelist$gene_symbol)]
DEG$Gene_symbol  = genelist$gene_symbol[match(rownames(DEG), genelist$ensembl_gene_id)]
write.csv(DEG, "./Data/DEG_Results_alpha05.csv")


# CNVs = makeGRangesFromDataFrame(top8)
# countOverlaps( LncRNA_Gr, CNVs)
# DEG$overlap[DEG$lncRNA] = countOverlaps( LncRNA_Gr, CNVs)

table(rownames(clinical)==rownames(datExpr0)) ## verify clinical matches expression for all samples
##  TRUE  
#    352 


gsg = goodSamplesGenes(datExpr0,verbose=4);
if (!gsg$allOK)
{
    if (sum(!gsg$goodGenes)>0)
        printFlush(paste("Removing genes:", paste(colnames(datExpr0)[!gsg$goodGenes], collapse = ", ")));
    if (sum(!gsg$goodSamples)>0)
        printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
    datExpr0= datExpr0[gsg$goodSamples, gsg$goodGenes]
    genelist=genelist[gsg$goodGenes,]
}

rm(Bspan_rows,SFARI,scores,Expr)


normadj <- (0.5+0.5*bicor(t(datExpr0)))^2
sdout <- 5
## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- ku-(mean(ku))/sqrt(var(ku))
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep=""))
print(colnames(t(datExpr0))[outliers])
print(table(outliers))

#"There are 0 outliers samples based on a bicor distance sample network connectivity standard deviation above 5"


datExpr0 <- log2(datExpr0+1)

powers = c(seq(8,14,by=1), seq(14,26, by=2));

##Save data for network construction on cluster
save.image("./Data/Pre_network_data.RData")

# #############################################################
####--  WARNING: The code within here requires a  large amount of RAM ( ~32GB)
#### we performed it on a linux machine with 16 cores ( if using windows, delete the nThreads arg in blockwiseModules)

# library(WGCNA)
# options(stringsAsFactors=FALSE)
# enableWGCNAThreads()
# 
# setwd("/scratch1/bgudena/WGCNA/ASD/")
# 
# 
# load(file="./Pre_network_data.RData")
# 
# sft=pickSoftThreshold(datExpr0, powerVector=powers, verbose=5, networkType="signed", corFnc = "bicor",corOptions = list(use = 'p', maxPOutliers = 0.1), blockSize = 30000)       
# 
# # sft=pickSoftThreshold(datExpr0, powerVector=powers, verbose=5, networkType="signed", corFnc = "bicor",corOptions = list(use = 'p', maxPOutliers = 0.1), blockSize = 30000)
#sft
# 
# net = blockwiseModules(datExpr0, power = 12,
#                        networkType="signed", minModuleSize = 50, maxBlockSize = 30000,
#                        mergeCutHeight = 0.15, deepsplit=4, corType= "bicor",corOptions = list(use = 'p', maxPOutliers = 0.1),
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "./TOM/TOM_ds4_min50",
#                        verbose = 3, nThreads = 16 )
# 
# table(net$colors)
# save.image("./Post_network_ds4_min50.RData")
##############################################################################

load(file="./Data/Post_network_ds4_min50.RData")

table(net$colors)

lncrnaColors = rep("grey", nrow(genelist))
DEGcolors = rep("grey", nrow(genelist))
ASDcolors = rep("grey", nrow(genelist))
lncrnaColors[genelist$lncRNA == TRUE] = "red3"
DEGcolors[genelist$L2FC != 0 ] = "blue3"
ASDcolors[!is.na(genelist$ASD_score) ] = "black"
moduleColors = labels2colors(net$colors)
genelist$Module = moduleColors

pdf("./Figures/Dendrogram.pdf")
plotDendroAndColors(net$dendrograms[[1]], cbind(moduleColors,lncrnaColors, DEGcolors, ASDcolors),
                    c("Module","LncRNAs","DEG", "ASD"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors

MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "./Data/networkconstruction.RData")

nGenes = ncol(datExpr0)
nSamples=nrow(datExpr0)


## calculate MEs with module colors
MEs0 <-moduleEigengenes(datExpr0,moduleColors)$eigengenes
MEs = orderMEs(MEs0)


counts = model.matrix( ~structure_acronym - 1, data = clinical)
colnames(counts) = str_replace(colnames(counts), pattern = "structure_acronym", replacement = "")


clinical = cbind(clinical, counts)
clinical$Male = rep(0, nrow(clinical))
clinical$Male[clinical$gender=="M"] = 1
    
    
moduleTraitCor = cor(MEs, clinical[ ,c(8:ncol(clinical))], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf("./Figures/Module-trait_relationships.pdf",width=14,height=7)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical[, c(8:ncol(clinical))]),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(60),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships") )

dev.off()

age=as.data.frame(clinical$Months.Post.Conception);
names(age)="Months_Age"
modNames=substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));


genes<-colnames(datExpr0)
genes2annot <- match(genes, genelist$ensembl_gene_id)
sum(is.na(genes2annot))


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, age, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(age), sep="");
names(GSPvalue) = paste("p.GS.", names(age), sep="");

geneInfo= data.frame(Ensembl_ID = genelist$ensembl_gene_id,
                     gene_symbol= genelist$gene_symbol,
                     Log2FC_DE = genelist$L2FC,
                     lncRNA=genelist$lncRNA,
                     ASD=genelist$ASD_score,
                     moduleColor = moduleColors,
                     geneTraitSignificance,
                     GSPvalue)
modOrder = order(-abs(cor(MEs, age, use = "p")));

for (mod in 1:ncol(geneModuleMembership))
{
    oldNames = names(geneInfo)
    geneInfo = data.frame(geneInfo, geneModuleMembership[, modOrder[mod]],
                          MMPvalue[, modOrder[mod]]);
    names(geneInfo) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                        paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
geneOrder = order(geneInfo$moduleColor, -abs(geneInfo$GS.Months_Age));
geneInfo = geneInfo[geneOrder, ]
write.csv(geneInfo, file = "./Data/geneInfo.csv",row.names = FALSE )

genelist$ASD_score = as.character(genelist$ASD_score)
genelist$ASD_score[is.na(genelist$ASD_score)]=0





library(dplyr)
moduleMembership = vector(mode = "numeric", length = nrow(geneInfo))
for (i in 1:nrow(geneInfo)){
    module = paste0("MM.",geneInfo$moduleColor[i])
    moduleMembership[i] = select(geneInfo,  match(module, colnames(geneInfo)))[i,]
}
geneInfo$Membership = moduleMembership

genelist$membership =geneInfo$Membership[match(genelist$ensembl_gene_id, geneInfo$Ensembl_ID)]


save.image("./Data/Post_geneinfo_network.RData")


##############################
library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("dark green","green","white","white","red","dark red"))(n = 599)


Module_dist = cor(MEs0,MEs0)
Modules0= paste(paste0(toupper(substr(names(table(moduleColors)), 1, 1)), tolower(substring(names(table(moduleColors)), 2))))
rownames(Module_dist) = Modules0
colnames(Module_dist) = Modules0

pdf("./Figures/ModuleCor_heatmap.pdf")
heatmap.2(Module_dist,trace="none",main="Module Correlation Matrix",RowSideColors=names(table(moduleColors)),notecol="black",key=TRUE,col=my_palette,symm=F,symkey=F,symbreaks=F)
dev.off()


###### this code is for modular analysis of ASD genes and lncRNAS
genelist$Module = moduleColors

genelist$ASD_score = as.character(genelist$ASD_score)
genelist$ASD_score[is.na(genelist$ASD_score)]=0

mod_sum<-matrix(data=NA,nrow=length(table(moduleColors)),ncol=4)
colnames(mod_sum) <- c("DE_lncRNAs","SFARI ASD Genes","ME16","total")
rownames(mod_sum) <- names(table(moduleColors))

ME16 = read.csv("./RAW/parikshak_m16.csv")
ME16 = ME16$ENSEMBL.GENE.ID
ME16_match = match(genelist$ensembl_gene_id, ME16)
genelist$ME16=0
genelist$ME16[!is.na(ME16_match)]=1
# LncRNAs <-as.list(table(factor(moduleColors[genelist$lncRNA== TRUE & genelist$L2FC == 0],lev=rownames(mod_sum))))
# DE_PC = table(factor(moduleColors[genelist$lncRNA== FALSE & genelist$L2FC != 0],lev=rownames(mod_sum)))
# mod_sum[,1] = as.numeric(DE_PC)

LncRNAs_DE = as.list(table(factor(moduleColors[genelist$lncRNA== TRUE & genelist$L2FC != 0],lev=rownames(mod_sum))))
mod_sum[,1] = as.numeric(LncRNAs_DE)


SFARI = table(factor(moduleColors[genelist$ASD_score=="1S" | genelist$ASD_score=="2S" | genelist$ASD_score=="3S" | genelist$ASD_score=="4S" | genelist$ASD_score==1 | genelist$ASD_score==2  | genelist$ASD_score==3 | genelist$ASD_score==4 | genelist$ASD_score==5],lev=rownames(mod_sum)))
mod_sum[,2] = as.numeric(SFARI)

#ASD_syndromic = as.list(table(factor(moduleColors[!is.na( genelist$ASD_score[ genelist$ASD_score=="1S" | genelist$ASD_score=="2S" | genelist$ASD_score=="3S" | genelist$ASD_score=="4S" | genelist$ASD_score=="S" ])], lev=rownames(mod_sum))))
ME16 = table(factor(moduleColors[genelist$ME16==1],lev=rownames(mod_sum)))
mod_sum[,3] = as.numeric(ME16)

mod_totals <- as.list(table(moduleColors))
mod_sum[,4]= as.numeric(mod_totals)

mod_sum = mod_sum[mod_sum[,1] > 1, ] ## remove modules with no lncRNAs


total_genes = nGenes
mat_p = matrix(data=NA, nrow=nrow(mod_sum), ncol=3)
mat_or = matrix(data=NA, nrow=nrow(mod_sum), ncol=3)
for (row in 1:nrow(mod_sum)){
    for (col in 1:3) {
        mod_total <- as.numeric(mod_sum[row, 4])
        mod_count <- as.numeric(mod_sum[row, col])
        mod_non <- mod_total-mod_count
        non_count <- sum(as.numeric(mod_sum[,col])) - mod_count       
        non_non <- (total_genes-mod_total)-non_count
        
        contigency <- matrix(c(mod_count, mod_non, non_count, non_non),2,2)
        colnames(contigency) = c(rownames(mod_sum)[row], paste0("Non-", rownames(mod_sum)[row]))
        rownames(contigency) = c(colnames(mod_sum)[col], paste0("Non-", colnames(mod_sum)[col]))
        
        
        results <- fisher.test(contigency, alternative = "greater")
        p_val <- results[[1]]
        OR <- results[[3]]
        
        mat_p[row,col]<-p_val
        mat_or[row,col]<-OR
    }
}

rownames(mat_p)<-rownames(mod_sum)
rownames(mat_or)<-rownames(mod_sum)
colnames(mat_p)<-colnames(mod_sum)[1:3]
colnames(mat_or)<-colnames(mod_sum)[1:3]


### now correct P-values for each test using FDR method
 ad_p = p.adjust(mat_p, method = "fdr")
dim(ad_p) = dim(mat_p)
colnames(ad_p) <-colnames(mat_p)
rownames(ad_p) <- rownames(mat_p)

#ad_p[ad_p < 10E-20] = 10E-20  ### adjust maximum bound of P-values to avoid saturation issues in heatmap

     ### -log10 transformation for figure

OR_filter<-matrix(data=NA,ncol=3,nrow=nrow(ad_p))
### * = p-value < 0.05 || ** = FDR-adjusted p-value < 0.05 ##### both need OR >= 1
OR_filter[ad_p <= 0.05] <- "*"
OR_filter[mat_or >= 1 & ad_p <= 0.05] <- paste0(round(mat_or[mat_or >= 1 & ad_p <= 0.05 ],2))

# log transform
ad_p<-log10(ad_p)*-1 

rownames(OR_filter)<-rownames(ad_p)



my_palette <- colorRampPalette(c("white","orange","red"))(n = 399)
 Modules_cap = paste0(toupper(substr(rownames(ad_p), 1, 1)), tolower(substring(rownames(ad_p), 2)))
rownames(ad_p) = Modules_cap

pdf("./Figures/module_geneset_enrichment.pdf")
heatmap.2(ad_p, cellnote=OR_filter, trace="none", main="ASD Gene Set Enrichment by Module", RowSideColors=Modules_cap, notecol="black",key=TRUE, col=my_palette, symm=F, symkey=F, symbreaks=F, key.xlab="-Log( FDR adjusted p-value )", notecex =1.5, margins =c(13,8), breaks = seq(0, 15, length.out =400))

dev.off()


# co-expression permutation testing ---------------------------------------

lncRNAlist = genelist[genelist$lncRNA & genelist$L2FC !=0, ]
#ASD_lncRNAs = lncRNAlist[lncRNAlist$Module=="blue" | lncRNAlist$Module=="magenta" | lncRNAlist$Module=="brown", ]



ASD_lncMatch = match( lncRNAlist$ensembl_gene_id , colnames(datExpr0))
genelist$ASD_score[is.na(genelist$ASD_score)] = 0
ASD_genes = genelist$ASD_score != 0 & genelist$ASD_score != "S" & genelist$ASD_score != "6" & genelist$ASD_score != "5" 

LncRNA_ASD_mat = bicor(datExpr0[ ,ASD_lncMatch], datExpr0[ ,ASD_genes], use="p", maxPOutliers = 0.1)

lncRNA_ASD_pairs=sum(abs(LncRNA_ASD_mat))


library(WGCNA)

P=10000 ## iterations
sig_pairs=c()
lnc_num = nrow(LncRNA_ASD_mat)
for (iter in 1:P){
    chosen = sample(1:ncol(datExpr0), lnc_num, replace = FALSE)
    chosen_ASD_mat = bicor( datExpr0[,chosen], datExpr0[ , ASD_genes], use="p", maxPOutliers = 0.1)
    sig_pairs=c(sig_pairs, sum(abs(chosen_ASD_mat)))
}



Z_scores = scale(c(sig_pairs,lncRNA_ASD_pairs)) ## calculate Z_scores of randomized distribution with actual value on end
Lnc_score = Z_scores[length(Z_scores)] ## Z_score of Actual lncRNA ASD pairs
pnorm(-abs(Lnc_score))

pdf("./Figures/LncRNA_coexpression_ASD-risk-genes.pdf")
hist(sig_pairs,,ylim=range(0,2600), xlab = "Summed Correlation", main = "LncRNAs co-expression to random gene-sets compared to ASD risk genes")
abline(v=lncRNA_ASD_pairs, lwd = 2, col = "red")
text(lncRNA_ASD_pairs-100,2500, pnorm(-abs(Lnc_score)), srt = 0.1, pos = 1, cex=1.2)
dev.off()


####################
#Perform same permutation test for ME16 genes
ME16 = read.csv("./RAW/parikshak_m16.csv")
ME16 = ME16$ENSEMBL.GENE.ID
ME16_match = match(genelist$ensembl_gene_id, ME16)
genelist$ME16=0
genelist$ME16[!is.na(ME16_match)]=1

ME16_genes = genelist$ME16==1

LncRNA_ASD_mat = bicor(datExpr0[,ASD_lncMatch], datExpr0[, ME16_genes], use="p", maxPOutliers = 0.1)

lncRNA_ASD_pairs=sum(abs(LncRNA_ASD_mat))


library(WGCNA)

P=10000 ## iterations
sig_pairs=c()
lnc_num = nrow(LncRNA_ASD_mat)
for (iter in 1:P){
    chosen = sample(1:ncol(datExpr0), lnc_num, replace = FALSE)
    chosen_ASD_mat = bicor( datExpr0[,chosen], datExpr0[,ME16_genes], use="p", maxPOutliers = 0.1)
    sig_pairs=c(sig_pairs, sum(abs(chosen_ASD_mat)))
}


Z_scores = scale(c(sig_pairs,lncRNA_ASD_pairs)) ## calculate Z_scores of randomized distribution with actual value on end
Lnc_score = Z_scores[length(Z_scores)] ## Z_score of Actual lncRNA ASD pairs
pval = pnorm(-abs(Lnc_score))


pdf("./Figures/LncRNA_coexpression_ME16-risk-genes.pdf")
hist(sig_pairs, xlim = range(14000,19000), xlab = "Summed Correlation", main = "LncRNAs co-expression to random gene-sets compared to ME16 genes")
abline(v= lncRNA_ASD_pairs, lwd = 2, col = "red")
text(lncRNA_ASD_pairs-100,1500, pval, srt = 0.1, pos = 1, cex=1.2)
dev.off()


DEG_match = match(rownames(DEG), genelist$ensembl_gene_id[ASD_genes])
DE_ASD = table(is.na(DEG_match))[1]
NotDE_ASD = length(genelist$ensembl_gene_id[ASD_genes]) - DE_ASD
DE_nonASD = nrow(DEG)-length(genelist$ensembl_gene_id[ASD_genes])
NotDE_nonASD = nrow(genelist) - nrow(DEG)

contingency = matrix(nrow=2,ncol=2, data = c(DE_ASD, DE_nonASD, NotDE_ASD, NotDE_nonASD))
rownames(contingency) = c("ASD","non-ASD")
colnames(contingency) = c("DE", "Not DE")
contingency
#         DE Not DE
# ASD       54    267
# non-ASD 1664  21035
fisher.test(contingency, "greater")
#  DEG enrichment of SFARI ASD risk genes
# Fisher's Exact Test for Count Data
# # 
# # data:  contingency
# # p-value = 1.344e-08
# # alternative hypothesis: true odds ratio is
# not equal to 1
# # 95 percent confidence interval:
# # 1.926684 3.622869
# # sample estimates:
# # odds ratio 
# # 2.664521 



load(file="./WGCNA/Data/Post_geneinfo_network.RData")
rm(datExpr0, MMPvalue, geneTraitSignificance, net)

#GTEx from 'http://www.gtexportal.org/static/datasets/gtex_analysis_v6/rna_seq_data/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.gct.gz'
GTEx = read.csv("./RAW/GTEx/GTEx_Analysis_v6_RNA-seq_RNA-SeQCv1.1.8_gene_median_rpkm.csv")
# GTEx_samples from http://www.gtexportal.org/home/tissueSummaryPage
GTEx_samples = read.csv("./RAW/GTEx//GTEx_sample_metadata.csv")
GTEx_samples = GTEx_samples[order(GTEx_samples$Tissue), ]  ## reorder samples to match GTEx columns (alphabetical)



match_GTEx = match(GTEx$Gene.Name, genelist$gene_symbol)
GTEx$module = genelist$Module[match_GTEx]
GTEx$lncRNA = genelist$lncRNA[match_GTEx]
GTEx$L2FC = genelist$L2FC[match_GTEx]

## remove genes not in network
GTEx = GTEx[!is.na(GTEx$module), ]


library(dplyr)
library(ggplot2)
library(stringr)
sub_GTEx = GTEx %>%
    group_by(module) %>%
    summarise(
        means = mean(Brain...Cortex)
    ) 
    ggplot(data = sub_GTEx, aes(module, means), color = module)+
    geom_bar( width=1, stat = "identity", aes(fill=module)) +
    scale_fill_identity()




exprMat = data.matrix(GTEx[, -c(1,2,56,57,58)])
rownames(exprMat) = GTEx$Gene.Name

##-- filter exprMat to remove tissues with less than 60 samples using GTEx_samples
exprMat = exprMat[GTEx$L2FC !=0 & GTEx$lncRNA , GTEx_samples$Number.of.RNASeq.Samples > 50]
# ### samples removed are
#     as.character(GTEx_samples$Tissue[ !GTEx_samples$Number.of.RNASeq.Samples > 50])
#     [1] "Minor Salivary Gland" "Kidney - Cortex"      "Bladder"              "Cervix - Ectocervix" 
#     [5] "Fallopian Tube"       "Cervix - Endocervix" 


brain_colors = rep("white", ncol(exprMat))
brain_colors[grepl("brain", str_to_lower(colnames(exprMat)))] = "blue"


exprMat = scale(t(exprMat))
exprMat = exprMat[ ,colSums(is.na(exprMat))!=nrow(exprMat)]

library(gplots)
library(RColorBrewer)
my_palette <- colorRampPalette(c("grey","white","orange","orangered","darkred","black"))(n = 599)
#DE_lncrns = scale(exprMat[GTEx$L2FC !=0 & GTEx$lncRNA, ])

pdf("./Figures/DE_lncRNAs_GTEX_heatmap.pdf", width = 16, height = 12)
#heatmap.2(exprMat[GTEx$L2FC !=0 & GTEx$lncRNA, ] , dendrogram = "column", trace="none", main="LncRNA Expression by Tissue-Type", notecol="black",key=TRUE, col=my_palette, symm=F, symkey=F, symbreaks=F, key.xlab="Scaled Median FPKM", scale ="row", margins = c(10,6))
heatmap.2(exprMat , dendrogram = "row", trace="none", main="LncRNA Expression by Tissue-Type", notecol="black", key=TRUE, col=my_palette, symm=F, symkey=F, symbreaks=F, key.xlab="Median FPKM Scaled by LncRNA",RowSideColors = brain_colors, margins = c(6,12))
dev.off()
#load(file="./Data/Post_geneinfo_network.RData")


load(file="./WGCNA/Data/Post_geneinfo_network.RData")

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("GOstats")
library(GOstats)
library(org.Hs.eg.db)
library(dplyr)


terms_full =data.frame()
for ( mod in unique(moduleColors)){
    
paramMF <- new("GOHyperGParams", geneIds = genelist$entrez_id[genelist$Module== mod], universeGeneIds = genelist$entrez_id, 
             ontology = c("MF"), annotation = "org.Hs.eg", pvalueCutoff = 0.01, testDirection = "over")

paramBP <- new("GOHyperGParams", geneIds = genelist$entrez_id[genelist$Module== mod], universeGeneIds = genelist$entrez_id, 
             ontology = c("BP"), annotation = "org.Hs.eg", pvalueCutoff = 0.01, testDirection = "over")
hypMF <- hyperGTest(paramMF)
hypBP <- hyperGTest(paramBP)

termsMF = summary(hypMF)
termsMF$type = "MF"
termsMF$Module =  mod
termsMF = rename(termsMF, GO_ID = GOMFID)

termsBP = summary(hypBP)
termsBP$type = "BP"
termsBP$Module =  mod
termsBP = rename(termsBP, GO_ID = GOBPID)

termsMF$Pvalue = p.adjust(termsMF$Pvalue, method = "fdr" , n = length(termsMF$Pvalue))
termsBP$Pvalue = p.adjust(termsBP$Pvalue, method = "fdr" , n = length(termsBP$Pvalue))

terms = rbind(termsMF[1:3, ], termsBP[1:3, ])

terms_full = rbind(terms_full, terms)
}


terms_full = terms_full[rowSums(is.na(terms_full)) != ncol(terms),  ]

write.csv(terms_full, "./WGCNA/Data/GOEnrichmentTable.csv")

terms_full = read.csv("./Data/GOEnrichmentTable.csv")

#terms_full = read.csv("./WGCNA/Data/GOEnrichmentTable.csv")
dir.create("./WGCNA/Figures/Expression_Trajectory")
library(ggplot2)
library(gridExtra)

for (mod in 1:ncol(MEs0)){
    modName = names(table(moduleColors))[mod]
    
    mod_terms = terms_full[terms_full$Module==modName, ]
    mod_terms = mod_terms[order(mod_terms$Pvalue), ]
    term_names = mod_terms$Term[order(mod_terms$Pvalue)]
    
    p1 = ggplot(mod_terms, aes( y=-log10(Pvalue), x=reorder(Term, Pvalue), fill= type)) +
        geom_bar( stat = "identity", position = position_dodge()) + theme_bw() +
        #theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        geom_hline(yintercept =  2, col = "red", lwd=1)+
        coord_flip()+
        ylab(expression(-Log[10]~Pvalue)) + xlab( "GO Term") +ggtitle(paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), "Functional Enrichment") )+ 
        scale_fill_manual(values=c("#56B4E9",  "#E69F00"))
    
    df = as.data.frame(cbind(Months=clinical$Months.Post.Conception[clinical$Months.Post.Conception <= 22], Expression = MEs0[,mod][clinical$Months.Post.Conception <= 22]) )
    p2= ggplot(df, aes(x = Months, y=Expression), main = paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), "Developmental Expression"),
               xlab="Post-Conception Developmental Months", ylab = "Module Eigengene (PC1)") +
        geom_point(shape=1, col ="grey") + theme_bw() +geom_jitter(alpha=0.7, shape=1)+
        geom_smooth(method=loess, size=1, col= modName, alpha=.5)+
        ggtitle(label = paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))) ) +
        geom_vline(xintercept =  10, col="blue", lwd=1)
    
    pdf(paste0("./WGCNA/Figures/Expression_Trajectory/", modName, ".pdf"), width = 12, height= 12)
    grid.arrange(p1, p2, ncol=2)
    
    dev.off()
}
