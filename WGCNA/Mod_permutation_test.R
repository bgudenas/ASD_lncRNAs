### Module characterization
    ## 5/17/16
    ## Brian Gudenas
    ## GOAL : Use randomized sampling to show modules have significant
    ##        Co-expression and DE
##_LIBRARY
library(dplyr)


workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)
load(file="./Data/Post_geneinfo_network.RData")
##--------------------------------------------------
genelist$Module = moduleColors

genesbyMod = group_by(genelist, Module)
Mod_DE = summarize(genesbyMod, Avg=mean(L2FC))
Mod_DE = Mod_DE[order(abs(Mod_DE$Avg), decreasing = TRUE),]


pdf("./Figures/DE_by_Module.pdf")
par(mar = c(6.5, 3, 3, 3));
df.bar = barplot(Mod_DE$Avg, col=Mod_DE$Module, ylab = "Log2 Fold Change", main="Differential expression of modules in ASD", las=2, names.arg = Mod_DE$Module)
lines(x= df.bar, y = as.numeric(rowMeans(Results)))
points(x= df.bar, y = as.numeric(rowMeans(Results)), col="red4", pch=16)
dev.off()
### Permutation testing



Iter = 10000

PermTest = function(Iter, genelist, test_var)
    #genesbyMod = group_by(genelist, Module)
    mod_counts = as.numeric(table(moduleColors))
    Results = data.frame(matrix(nrow=length(mod_counts), ncol = Iter, data=0))
        rownames(Results) = names(table(moduleColors))

    for (I in 1:Iter){
        avail_genes = 1:nrow(genelist)   ##establish available genes
        for (mod in 1:length(mod_counts)) {
            rand_genes = sample(avail_genes, mod_counts[mod])
                avail_genes=avail_genes[-(rand_genes)]   ## remove chosen genes from available
            rand_avg = mean(genelist$L2FC[rand_genes])
            Results[mod,I] = rand_avg
            
                
        }
    }
    
### statistical Enrichment?
Z_scores = scale(c(as.numeric(Results[2,]),Mod_DE$Avg[2]))


#Co-expression
##########################################
Real_coexp = data.frame(matrix(nrow=length(mod_counts), ncol = 1, data=0))
rownames(Real_coexp) = names(table(genelist$Module))
for (mod in names(table(genelist$Module))) {

    colnames(Real_coexp) = "Mean Module Correlation"
    actual_coexp = mean(cor(datExpr0[,genelist$Module == mod], datExpr0[,genelist$Module == mod]))
    Real_coexp[rownames(Real_coexp)==mod,1] = actual_coexp
}
    




cor_mat=cor(datExpr0,datExpr0)
rm(GSPvalue, datExpr0, geneModuleMembership, geneTraitSignificance, net, lnc_test, geneTree, MMPvalue, MEs, MEs0)


permTest = function() {
Iter = 1000
mod_counts = as.numeric(table(moduleColors))
Results = matrix(nrow=length(mod_counts), ncol = Iter, data=0)
rownames(Results) = names(table(moduleColors))

for (I in 1:Iter){
    avail_genes = 1:length(moduleColors)
    for (mod in 1:length(mod_counts)) {
        rand_genes = sample(avail_genes, mod_counts[mod])
        avail_genes = avail_genes[-rand_genes]   ## remove chosen genes from available
        
        rand_sum = .Internal(mean(cor_mat[rand_genes,rand_genes]))
        #rand_sum = mean(cor(datExpr0[,rand_genes], datExpr0[,rand_genes]))
        Results[mod,I] = rand_sum 
    }
}
#return(Results)
#}
#system.time(permTest())
#1:24 @ 5k    

rm(cor_mat)
save.image("./Data/CoEXP_perm_results.RData")
