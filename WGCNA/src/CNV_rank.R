load(file="./WGCNA/Data/Post_geneinfo_network.RData")
library(GenomicRanges)
library(WGCNA)
CNVs = read.csv("./WGCNA/RAW/SFARI/cnv-summary.csv")
CNVs = CNVs[CNVs$Report.Class == "Major", ]

source("./cyto_converter/cyto_converter.R")
CNV_df = cyto_converter(CNVs$CNV.Locus, cytobands = cytobands)
table(rowSums(is.na(CNV_df)))

CNV_df = CNV_df[rowSums(is.na(CNV_df))==0, ]


CNVs = makeGRangesFromDataFrame(CNV_df)

genelist_Gr = makeGRangesFromDataFrame(genelist[!is.na(genelist$start), colnames(genelist)=="chromosome" | colnames(genelist)=="start" |  colnames(genelist)=="end" ], start.field = "start", end.field = "end" )
names(genelist_Gr) = genelist$gene_symbol[!is.na(genelist$start) ]
countOverlaps( LncRNA_Gr, CNVs)

#genelist_Gr = makeGRangesFromDataFrame(genelist[ !is.na(genelist$start), c(13, 14, 15)], start.field = "start", end.field = "end" )
 
#  DEG$overlap[DEG$lncRNA] = countOverlaps( LncRNA_Gr, CNVs)
# DEG = DEG[order(DEG$overlap, decreasing = TRUE), ]


overlaps = countOverlaps( genelist_Gr, CNVs)

genelist$overlaps = overlaps[match(genelist$gene_symbol, names(overlaps))]
#genelist$overlaps[is.na(genelist$overlaps)] = 0
library(dplyr)
library(ggplot2)

mod_overlaps = genelist %>%
    group_by(Module) %>% 
    summarise(avg_overlap = mean(overlaps)) %>% 
    arrange(desc(avg_overlap) ) %>% 
    ggplot()+
    geom_bar(aes(x=Module, y = avg_overlap, fill=Module), stat = "identity") +
    scale_fill_identity() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))


ranked = genelist[genelist$lncRNA == TRUE & genelist$L2FC !=0, ]
ranked = ranked[ranked$Module == "red" | ranked$Module == "blue" | ranked$Module == "brown", ] 

ranked$nearest = DEG$nearest[match(ranked$ensembl_gene_id, rownames(DEG))]
ranked$nearest_L2FC = genelist$L2FC[match(ranked$nearest, genelist$gene_symbol)]

#ranked$nearest_L2FC = DEG$log2FoldChange[match(DEG$nearest, rownames(DEG))]
ranked$nearest_ASD_score = DEG$nearest_ASD_score[match(ranked$ensembl_gene_id, rownames(DEG))]

  # write.csv(ranked, "./Data/ranked_lncRNAs.csv")
    
#     
# genelist %>% 
#     group_by(ASD_score) %>% 
#     summarise(avg = mean(overlaps))


#ranked = read.csv("./WGCNA/Data/Ranked_lncRNAs.csv")

top_cors = data.frame()
for ( i in 1:nrow(ranked)){
    lncRNA = ranked$ensembl_gene_id[i]
    
    lnc_expr = datExpr0[ ,colnames(datExpr0)==lncRNA]
    
    cor_mat = bicor(datExpr0[ ,colnames(datExpr0)!=lncRNA & genelist$biotype=="protein_coding"], lnc_expr,  use="p", maxPOutliers = 0.1)
    #which.max(cor_mat)
    ind = cor_mat[which.max(cor_mat), ]

    gene_name = genelist[genelist$ensembl_gene_id == names(ind),  ]
    E_ID = genelist$ensembl_gene_id[match(ranked$nearest[i], genelist$gene_symbol)]
    near_cor = bicor(datExpr0[ ,colnames(datExpr0) == E_ID], lnc_expr,  use="p", maxPOutliers = 0.1)
    gene_name$nearest_cor = near_cor
    
    gene_name$Top_cor = as.numeric(ind)
    top_cors = rbind(top_cors, gene_name)
}
ranked = cbind(ranked,top_cors)

ranked = ranked[order(ranked$overlaps, decreasing = TRUE), ]

write.csv(ranked, "./WGCNA/Data/ranked_lncRNAs.csv")

ranked = read.csv("./WGCNA/Data/Ranked_lncRNAs.csv")

contigency = matrix(nrow=2, ncol=2)
rownames(contigency) = ""
    