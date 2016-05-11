#5/3/16
#### load WGCNA v1.46
library(WGCNA)
options(stringsAsFactors=FALSE)
enableWGCNAThreads()

workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)


Bspan_rows<-read.csv("./RAW/Brainspan/rows_metadata.csv")  ##genelist of Brainspan transcriptome summarized to genes,5/15
DEG <-read.csv("C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/DESeq2/DEG_Results_alpha05.csv",row.names = 1) ## ASD Differentially expression 


SFARI <-read.csv("./RAW/SFARI/gene-summary.csv") #7/11/15 SFARI ASD genes, 707 total
scores<-read.csv("./RAW/SFARI/gene-score.csv") #SFARI ASD risk gene scores

score_match = match(scores$Gene.Symbol, SFARI$Gene.Symbol)
scores$entrez = SFARI$Entrez.GeneID[score_match]


## load expression matrix and sample metadata
Expr <- read.csv("./Raw/Brainspan/expression_matrix.csv",header=FALSE,row.names=1)  ##RPKM expression data
clinical <- read.csv("./Raw/Brainspan/Bspan_clinical.csv",row.names=1)  ##matrix containing sample metadata
Expr=Expr[ , clinical$Temporal==1 | clinical$Neocortex==1]
clinical = clinical[clinical$Temporal==1 | clinical$Neocortex==1,]
colnames(Expr)<-rownames(clinical)


#datExpr0<-t(Expr[unique(c(gene_list$row_num, ASD_nums)),])
#colnames(datExpr0) <-unique(c(as.character(gene_list$ensembl_gene_id), as.character(ASD_names)))
datExpr0 = t(Expr)
colnames(datExpr0)=Bspan_rows$ensembl_gene_id

#-- THis loop checks the # of times each gene is expressed in a sample >= 1 RPKM
vals=c()
for (i in 1:ncol(datExpr0)) {
    vals=c(vals,sum(datExpr0[,i] >= 1) )
}

Gene_filter = (vals > 3 )  ## genes required to have normalized RPKM >=1 in atleast 3/331 of samples
datExpr0 = datExpr0[, Gene_filter]


##--Create Genelist dataframe containing all gene info from (DEG and SFARI)
genelist = Bspan_rows[Gene_filter,]
DEG_map = match(genelist$ensembl_gene_id, rownames(DEG))


genelist$L2FC = DEG$log2FoldChange[DEG_map]
genelist$L2FC[is.na(genelist$L2FC)] = 0    ### Give 0 to any gene not Differentially expressed

mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","gene_biotype","chromosome_name","start_position","end_position"),
            filters = "ensembl_gene_id", values =  genelist$ensembl_gene_id)

mart_map = match(genelist$ensembl_gene_id, map$ensembl_gene_id)

genelist$biotype = map$gene_biotype[mart_map]

ASD_match = match(genelist$entrez_id, scores$entrez)
genelist$ASD_score = scores$Score[ASD_match]



table(rownames(clinical)==rownames(datExpr0)) ## verify clinical matches expression for all samples
##  TRUE  
#    331


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

rm(Bspan_rows,SFARI,scores,Expr,map)


datExpr0<-log2(datExpr0+1)

powers = c(seq(8,14,by=1), seq(14,26, by=2));

##Save data for network construction on cluster
save.image("Pre_network_data.RData")

# #############################################################
# ##-- all code within here is used in create_network.R for cluster use
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
# 
# net = blockwiseModules(datExpr0, power = 13,
#                        networkType="signed", minModuleSize = 20, maxBlockSize = 30000,       
#                        mergeCutHeight = 0.2,deepsplit=4, corType= "bicor",corOptions = list(use = 'p', maxPOutliers = 0.1),
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "./TOM/TOM",
#                        verbose = 3 )
# 
# table(net$colors)
# save.image("./Post_network.RData")
##############################################################################

load(file="./Data/Post_network.RData")

lncRNA_filter = c("3prime_overlapping_ncrna", "antisense", "lincRNA", "processed_transcript", "sense_intronic" , "sense_overlapping")

lnc_test=c()
for (biotype in genelist$biotype) {
    lnc_test =c(lnc_test, (sum(grepl(biotype, lncRNA_filter)) > 0))
}
genelist$lncRNA = lnc_test



table(net$colors)
pdf("./Figures/Dendrogram.pdf")
moduleColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], moduleColors,
                    c("Module","Type"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "./Data/networkconstruction.RData")
nGenes<-ncol(datExpr0)
nSamples=nrow(datExpr0)


## calculate MEs with module colors
MEs0 <-moduleEigengenes(datExpr0,moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, clinical[,c(1:2)], use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

pdf("./Figures/Module-trait_relationships.pdf",width=14,height=7)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));

labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(clinical[, c(1:2)]),
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

age=as.data.frame(clinical$MonthsAge);
names(age)="Months_Age"
modNames=substring(names(MEs),3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));


genes<-colnames(datExpr0)
genes2annot <- match(genes,genelist$ensembl_gene_id)
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

save.image("./Data/Post_geneinfo_network.RData")

load(file="./Data/Post_geneinfo_network.RData")


##############################
library(gplots)
library(RColorBrewer)
###### this code is for modular analysis of ASD genes and lncRNAS
genelist$ASD_score = as.character(genelist$ASD_score)
genelist$ASD_score[is.na(genelist$ASD_score)]=0

mod_sum<-matrix(data=NA,nrow=length(table(moduleColors)),ncol=4)
colnames(mod_sum) <- c("LncRNAs_DE","ASD_high","ASD_Full","total")
rownames(mod_sum) <- as.list(names(table(moduleColors)))

# LncRNAs <-as.list(table(factor(moduleColors[genelist$lncRNA== TRUE & genelist$L2FC == 0],lev=rownames(mod_sum))))
# DE_PC = table(factor(moduleColors[genelist$lncRNA== FALSE & genelist$L2FC != 0],lev=rownames(mod_sum)))
# mod_sum[,1] = as.numeric(DE_PC)

LncRNAs_DE = as.list(table(factor(moduleColors[genelist$lncRNA== TRUE & genelist$L2FC != 0],lev=rownames(mod_sum))))
mod_sum[,1] = as.numeric(LncRNAs_DE)

ASD_high = table(factor(moduleColors[genelist$ASD_score==1 | genelist$ASD_score==2  | genelist$ASD_score==3 | genelist$ASD_score==4 ],lev=rownames(mod_sum)))
mod_sum[,2] = as.numeric(ASD_high)

#ASD_syndromic = as.list(table(factor(moduleColors[!is.na( genelist$ASD_score[ genelist$ASD_score=="1S" | genelist$ASD_score=="2S" | genelist$ASD_score=="3S" | genelist$ASD_score=="4S" | genelist$ASD_score=="S" ])], lev=rownames(mod_sum))))
ASD_Full = table(factor(moduleColors[genelist$ASD_score=="1S" | genelist$ASD_score=="2S" | genelist$ASD_score=="3S" | genelist$ASD_score=="4S" | genelist$ASD_score==1 | genelist$ASD_score==2  | genelist$ASD_score==3 | genelist$ASD_score==4 | genelist$ASD_score == 5 ],lev=rownames(mod_sum)))
mod_sum[,3] = as.numeric(ASD_Full)

mod_totals <- as.list(table(moduleColors))
mod_sum[,4]= as.numeric(mod_totals)

mod_sum = mod_sum[mod_sum[,1] > 0,] ### this removes all modules without lncRNAs to lower the # of statiscal tests sense we are only testing for alternative = greater

total_genes = nGenes
mat_p = matrix(data=NA, nrow=nrow(mod_sum), ncol=3)
mat_or = matrix(data=NA, nrow=nrow(mod_sum), ncol=3)
for (row in 1:nrow(mod_sum)){
    for (col in 1:3) {
        mod_total <- as.numeric(mod_sum[row,4])
        mod_count <- as.numeric(mod_sum[row,col])
        mod_non <- mod_total-mod_count
        non_count <- sum(as.numeric(mod_sum[,col])) - mod_count       
        non_non <- (total_genes-mod_total)-non_count
        
        contigency <- matrix(c(mod_count, mod_non, non_count, non_non),2,2)
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

ad_p<-c()

### now correct P-values for each test using FDR method
for (i in 1:ncol(mat_p)){
    ad_p<-cbind(ad_p,p.adjust(mat_p[,i],"fdr"))
}
colnames(ad_p) <-colnames(mat_p)

log_p<-log10(mat_p)*-1      ### -log10 transformation for figure

my_palette <- colorRampPalette(c("white","orange","red"))(n = 299)
col_breaks = c(seq(0,1.2,length=100), # for red
               seq(1.2,2,length=100), # for yellow
               seq(2,6,length=100)) # for green

OR_filter<-matrix(data=NA,ncol=3,nrow=nrow(log_p))
OR_filter[mat_or >= 1 &log_p >= 1.3]<-round(mat_or[mat_or >= 1 & log_p>= 1.3],2)

rownames(OR_filter)<-rownames(log_p)

pdf("./Figures/concise_module_geneset_enrichment.pdf")
par(mar=c(5,2,3,5))

heatmap.2(log_p,cellnote=OR_filter,trace="none",main="Risk gene Enrichment",RowSideColors=rownames(log_p),notecol="black",key=TRUE,col=my_palette,symm=F,symkey=F,symbreaks=F)
#heatmap.2(log_p[rowSums(OR_filter[,-c(1,6)],na.rm=TRUE)!=0,-c(1,6)],cellnote=OR_filter[rowSums(OR_filter[,-c(1,6)],na.rm=TRUE)!=0,-c(1,6)],margins=c(8,6),trace="none",main="Risk gene Enrichment",RowSideColors=rownames(log_p[rowSums(OR_filter[,-c(1,6)],na.rm=TRUE)!=0,-c(1,6)]),notecol="black",key=TRUE,col=my_palette,symm=F,symkey=F,symbreaks=F)
dev.off()



###9/9 CNV detection
setwd("ASD")
load(file="./Data/CNVcleaning.RData")
rm(GSPvalue,OR_filter,OR,ad_p,cont_CNVs,geneModuleMembership,mod_sum,mod_count,mod_totals,textMatrix,ensembl_hg36,col_breaks,p_val,sft,powers,my_palette,MMPvalue,gene_matches,mat_p,mat_or,geneTraitSignificance,gsg,mod,gene,net,results,geneTree,mod_non)

rawCNVs <-read.csv("Raw/SFARI/individual-data.csv") ##7/24
SFARI_cyto <- read.csv("CNV/SFARI_converted_hg38.csv") ###cytoband converted to hg38 coordinates using cyto_converter.R
rawCNVs=rawCNVs[complete.cases(SFARI_cyto),]
SFARI_cyto=SFARI_cyto[complete.cases(SFARI_cyto),]

SFARI_cyto$chromosome =paste("chr",SFARI_cyto$chromosome,sep="")

library(GenomicRanges)
G_SFARI=makeGRangesFromDataFrame(SFARI_cyto[,2:4])
names(G_SFARI)= SFARI_cyto$cytoband

ASD_CNVs = G_SFARI[rawCNVs$Primary.Diagnosis=="Autism" | rawCNVs$Primary.Diagnosis=="ASD" | rawCNVs$Primary.Diagnosis=="Aspergers",]
cont_CNVs = G_SFARI[rawCNVs$Primary.Diagnosis=="Control" | rawCNVs$Primary.Diagnosis=="Control (matched sibling)" ,]

save.image("./Data/CNVcleaning.RData")


