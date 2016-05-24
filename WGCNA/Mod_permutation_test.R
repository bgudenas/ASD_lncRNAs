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
#Mod_DE = Mod_DE[order(abs(Mod_DE$Avg), decreasing = TRUE),]

### Permutation testing
Iter = 10000
    #genesbyMod = group_by(genelist, Module)
    mod_counts = as.numeric(table(moduleColors))
    Results = matrix(nrow=length(mod_counts), ncol = Iter, data=0)
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

Pvals =c()
for (i in 1:nrow(Results)){
    ### statistical Enrichment
    Z_scores = scale(c(as.numeric(Results[i, ]), as.numeric(Mod_DE[i,2]) ))   ## transform all permuted and Real data to Z-scores
    Real = tail(Z_scores, 1)
    Pvals=c(Pvals, 2*pnorm(-abs(Real)))   ## use pnorm to calculate two-sided P-value of real Z-score and append
}
Pvals = Pvals*length(Pvals) ## adjust P-values for multiple testing

Sig_P = rep("", length(Pvals))
Sig_P[Pvals < 0.05]="*"


## reorder Sig_P, Mod_DE and Results by Mod_DE
Results = Results[order(abs(Mod_DE$Avg), decreasing = TRUE),]
Sig_P = Sig_P[order(abs(Mod_DE$Avg), decreasing = TRUE)] 
Mod_DE = Mod_DE[order(abs(Mod_DE$Avg), decreasing = TRUE),] 

Sig_pos = rep(1, length(Sig_P))
Sig_pos[Mod_DE$Avg > 0] =3

pdf("./Figures/DE_by_Module.pdf")
par(mar = c(6.5, 4, 3, 3));
df.bar = barplot(Mod_DE$Avg, col=Mod_DE$Module, ylab = "Log2 Fold Change", main="Differential expression of modules in ASD", las=2, names.arg = Mod_DE$Module)
lines(x= df.bar, y = as.numeric(rowMeans(Results)))
points(x= df.bar, y = as.numeric(rowMeans(Results)), col="red4", pch=16)
text(x= df.bar, 0, Sig_P, cex =2, pos=Sig_pos, col="white")
## white asterisks denote significance at p = 0.05
dev.off()


#Co-expression
##########################################
Real_coexp = data.frame(matrix(nrow=length(mod_counts), ncol = 1, data=0))
rownames(Real_coexp) = names(table(genelist$Module))
for (mod in names(table(genelist$Module))) {

    colnames(Real_coexp) = "Correlation"
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

workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)
load(file="Data/CoEXP_perm_results.RData")

Pvals =c()
for (i in 1:nrow(Results)){
    ### statistical Enrichment
    Z_scores = scale(c(as.numeric(Results[i,]),Real_coexp[i,]))   ## transform all permuted and Real data to Z-scores
    Real = tail(Z_scores, 1)
    Pvals=c(Pvals, 2*pnorm(-abs(Real)))   ## use pnorm to calculate P-value of real Z-score and append
   }
Pvals = Pvals*length(Pvals) ## adjust P-values for multiple testing
## all modules signficiant at P = .001 Except the Gray Module

## create vector to represent signficance as *
Sig_P = rep("", length(Pvals))
Sig_P[Pvals < 0.001]="*"    



pdf("./Figures/CoXP_perm.pdf")
par(mar = c(6.5, 5, 3, 3));
df.bar = barplot(Real_coexp$Correlation, col=rownames(Real_coexp), ylab = "Pearson Correlation Coefficient", main="Average Module Coexpression", las=2, names.arg = rownames(Real_coexp))
lines(x= df.bar, y = as.numeric(rowMeans(Results)), lwd=2)
points(x= df.bar, y = as.numeric(rowMeans(Results)), col="red4", pch=16)
#text(x= df.bar, 0, Sig_P, cex =2, pos=3, offset = 14, col="green4")
dev.off()
