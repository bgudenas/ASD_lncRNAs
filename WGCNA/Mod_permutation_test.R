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
    

Z_scores = scale(c(as.numeric(Results[2,]),Mod_DE$Avg[2]))
    