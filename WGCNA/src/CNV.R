# CNV analysis ------------------------------------------------------------
library(biomaRt)
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(broom)
library(WGCNA)
options(stringsAsFactors=FALSE)


workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)
load(file="./Data/Post_geneinfo_network.RData")
rm( MMPvalue, geneTraitSignificance, net, GSPvalue, geneModuleMembership, geneTree, MEs0, textMatrix, moduleTraitCor, moduleTraitPvalue, gsg, i, MEs)


mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
map = getBM(mart = mart, attributes = c("ensembl_gene_id","band","gene_biotype","chromosome_name","start_position","end_position"),
            filters = "ensembl_gene_id", values =  genelist$ensembl_gene_id)

mart_map = match(genelist$ensembl_gene_id, map$ensembl_gene_id)

genelist$Chromosome = map$chromosome_name[mart_map]
genelist$Start = map$start_position[mart_map]
genelist$End = map$end_position[mart_map]

genelist = genelist[rowSums(is.na(genelist[ ,c(11,12,13)])) == 0, ] ## remove 1618 genes without locations
genelist$Width = genelist$End - genelist$Start

CNVs = read.csv("RAW/SFARI/cnv-summary.csv") ## from SFARI with added column for row.nums
CNVs = CNVs[CNVs$Report.Class=="Major", ]

source("../cyto_converter/cyto_converter.R") ## functions to convert cytoband to hg38 genomic coords
CNVs_coords = cyto_converter(bands = CNVs$CNV.Locus, cytobands = cytobands)
CNVs_coords$Type = CNVs$CNV.Type
CNVs_coords$Report.Class = CNVs$Report.Class 
CNVs_coords = CNVs_coords[rowSums(is.na(CNVs_coords)) == 0, ] ## remove 5 cytobands unable to be converted


#  Make GRanges for each CNV type -----------------------------------------
gDups = makeGRangesFromDataFrame(CNVs_coords[CNVs_coords$Type == "Duplication", ])
gDels = makeGRangesFromDataFrame(CNVs_coords[CNVs_coords$Type == "Deletion", ])
gCNVs = makeGRangesFromDataFrame(CNVs_coords)

gGenes = makeGRangesFromDataFrame(genelist[ , c(11,12,13)] )
names(gGenes) = genelist$ensembl_gene_id

overlaps_dups = countOverlaps(query = gGenes, subject = gDups, type = "within")
overlaps_dels = countOverlaps(query = gGenes, subject = gDels, type = "within")
overlaps_CNVs = countOverlaps(query = gGenes, subject = gCNVs, type = "within")

genelist$Dups = overlaps_dups
genelist$Dels = overlaps_dels
genelist$CNVs = overlaps_CNVs




# attach module membership vals -------------------------------------------

moduleMembership = vector(mode = "numeric", length = nrow(geneInfo))
for (i in 1:nrow(geneInfo)){
    module = paste0("MM.",geneInfo$moduleColor[i])
    moduleMembership[i] = as.numeric(select_(geneInfo[i,], module))
}
geneInfo$Membership = moduleMembership

genelist$membership =geneInfo$Membership[match(genelist$ensembl_gene_id, geneInfo$Ensembl_ID)]


genelist$DE_lncRNA = genelist$L2FC!=0 & genelist$lncRNA
# table(genelist$DE_lncRNA)

# genelist = genelist[order(genelist$CNVs, decreasing = TRUE), ]



ranked = genelist[genelist$DE_lncRNA, ]

ASDranked = read.csv("../Manuscript/ASD_lncRNAs_coexpression_draft3/Table S3. Ranked_lncRNAs.csv")
ranked$ASD_sumCor = ASDranked$ASD_sumCor[match(ranked$ensembl_gene_id, ASDranked$ensembl_gene_id)]
ranked = ranked %>%
    arrange(desc(ASD_sumCor))


genelist$Class = "Background"
genelist$Class[genelist$L2FC!=0 & !genelist$lncRNA] = "DE_PC"
genelist$Class[genelist$ASD_score=="1S" | genelist$ASD_score=="2S" | genelist$ASD_score=="3S" | genelist$ASD_score=="4S" | genelist$ASD_score==1 | genelist$ASD_score==2  | genelist$ASD_score==3 | genelist$ASD_score==4 | genelist$ASD_score==5] = "SFARI"
genelist$Class[genelist$L2FC!=0 & genelist$lncRNA  ] = "DE_lncRNA"




tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "SFARI"], y = genelist$CNVs[genelist$Class == "Background"], paired = FALSE, alternative = "greater"))
# statistic      p.value                                            method alternative
# 1   4508532 4.384798e-05 Wilcoxon rank sum test with continuity correction     greater

tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "DE_lncRNA"], y = genelist$CNVs[genelist$Class == "Background"], paired = FALSE, alternative = "greater"))
tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "SFARI"], y = genelist$CNVs[genelist$Class == "DE_lncRNA"], paired = FALSE))



top5 = vector(mode="character", length=nrow(ranked))
datExpr0 = as.data.frame(datExpr0)
datExpr = datExpr0[ , !is.na(match(colnames(datExpr0), genelist$ensembl_gene_id))] ## remove genes without locations
for ( i in 1:nrow(ranked)) {
    E_id = as.character(ranked$ensembl_gene_id[i])
    corrs = bicor(select_(datExpr, E_id), datExpr)
    corrs = corrs[ ,order(corrs, decreasing = TRUE) ]
    top = names(corrs)[c(2:6)]
    
    top5[i] = list(genelist$gene_symbol[match(top, genelist$ensembl_gene_id)])
}

ranked = cbind(ranked, data.frame(matrix(unlist(top5), nrow=129, byrow=T)))
#colnames(ranked[,(ncol(ranked)-4):(ncol(ranked))]) = c("TopCor1","TopCor2","TopCor3", "TopCor4","TopCor5")
    
write.csv(ranked, "./Data/Ranked_lncRNAs.csv")

# Exploratory data analysis -----------------------------------------------



genelist %>%
    group_by(Class) %>%
    summarise(sumOverlap = sum(CNVs), averageOverlap =mean(CNVs), medianCNV = median(CNVs), L2FC = mean(L2FC, na.rm=TRUE), aveMembership = mean(membership), count =n()) %>%
    write.csv("./Data/CNVs_by_Class.csv")

CNVs_byMod = genelist %>%
    group_by(Module) %>%
    summarise(sumOverlap = sum(CNVs), averageOverlap = mean(CNVs), averageDups = mean(Dups), averageDels = mean(Dels), aveMembership = mean(membership), count =n() ) %>%
    arrange(desc(averageOverlap))
    



bg_genelist = genelist[, -20]
    ggplot(genelist, aes(x= membership, fill = Class), ) +
    geom_histogram(data = bg_genelist, fill= "grey", alpha = .5) +
    geom_histogram(color = "black") +
    facet_wrap(~ Class) +
    
    theme_bw()




genelist %>%
    group_by(Class) %>%
    ggplot( aes(y=CNVs, x = Class)) +
    geom_boxplot( aes(fill=Class)) 

genelist %>%
    group_by(Module) %>%
    #summarise(averageOverlap =mean(Overlaps)) %>%
    ggplot(aes(y=CNVs, x=Module)) +
    geom_boxplot( aes(fill=Module)) +
    scale_fill_identity() +
    geom_hline(yintercept = median(genelist$CNVs), col = "red")


wilc_ASD.nonASD = tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "SFARI"], y = genelist$CNVs[genelist$Class != "SFARI"], paired = FALSE, alternative = "greater"))
# statistic    p.value                                            method alternative
# 1   4669447 0.01159776 Wilcoxon rank sum test with continuity correction     greater


tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "SFARI"], y = genelist$CNVs[genelist$Class== "DE_lncRNA"], paired = FALSE, alternative = "greater"))
tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "SFARI"], y = genelist$CNVs[genelist$Class== "DE_PC"], paired = FALSE, alternative = "greater"))
tidy(wilcox.test(x = genelist$CNVs[genelist$Class == "DE_lncRNA"], y = genelist$CNVs[genelist$Class== "Background"], paired = FALSE, alternative = "greater"))
# statistic    p.value                                            method alternative
# 1   4664501 0.01285282 Wilcoxon rank sum test with continuity correction     greater
# ASD genes are within CNVs at a greater degree than non-ASD genes


tidy(wilcox.test(x = ranked$Dups, y = ranked$Dels, paired = TRUE))

tidy(wilcox.test(x = ranked$Dups, y = ranked$Dels, paired = TRUE))


ranked_summary = ranked %>%
    group_by(Module) %>%
    summarise(totalOverlap =sum(CNVs),meanOverlap = mean(CNVs), L2FC = mean(L2FC, na.rm=TRUE), aveMembership = mean(membership), count= n(), width = sum(Width) ) %>%
    arrange(desc(totalOverlap))

bg_lncRNAs = genelist[genelist$lncRNA & genelist$L2FC == 0, ] %>%
    group_by(Module) %>%
    summarise(totalOverlap =sum(CNVs), meanOverlap = mean(CNVs), L2FC = mean(L2FC, na.rm=TRUE), aveMembership = mean(membership), count= n(), width = sum(Width) ) %>%
    arrange(desc(totalOverlap))
bg_lncRNAs$aveSize =  bg_lncRNAs$width / bg_lncRNAs$count




ranked %>%
    group_by(biotype) %>%
    summarise( averageOverlap =mean(CNVs), medianCNV = median(CNVs), L2FC = mean(L2FC, na.rm=TRUE), aveMembership = mean(membership), count = n()) %>%
    arrange(desc(averageOverlap))


