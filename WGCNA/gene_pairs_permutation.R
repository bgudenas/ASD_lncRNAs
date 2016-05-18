##-- ASD-genepair permutation test
# Ho = DE_lncRNAs have average coexpression with ASD risk genes
# HA = DE_lncRNAs have greater coexpression with ASD risk genes than average

workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)

load(file="./Data/Post_geneinfo_network.RData")

DE_lncRNAs = genelist$lncRNA & genelist$L2FC !=0 ## must be lncRNA and DE
ASD_genes = genelist$ASD_score != 0 & genelist$ASD_score != "S" & genelist$ASD_score != "6" & genelist$ASD_score != "5"

LncRNA_ASD_mat = cor(datExpr0[,DE_lncRNAs], datExpr0[,DE_lncRNAs])

lncRNA_ASD_pairs=sum(LncRNA_ASD_mat > 0.8)

### permutation loop
P=10000 ## iterations

#fExpr = datExpr0[,-c(DE_lncRNAs, ASD_genes) ]  ## eligible genes
sig_pairs=c()
lnc_num = sum(DE_lncRNAs)
for (iter in 1:P){
    chosen = sample(1:ncol(datExpr0), lnc_num)
    chosen_ASD_mat = cor( datExpr0[,chosen], datExpr0[,chosen])
    sig_pairs=c(sig_pairs, sum(chosen_ASD_mat > 0.8))
}

Z_scores = scale(c(sig_pairs,lncRNA_ASD_pairs)) ## calculate Z_scores of randomized distribution with actual value on end
Lnc_score = Z_scores[length(Z_scores)] ## Z_score of Actual lncRNA ASD pairs
pnorm(-abs(Lnc_score))
#p value = 1.618207e-26

pop_sd =  sd(c(sig_pairs,lncRNA_ASD_pairs))
pop_mean = mean(c(sig_pairs,lncRNA_ASD_pairs))
z = (lncRNA_ASD_pairs - pop_mean) / (pop_sd/ sqrt(P))
pnorm(-abs(z))


pdf("./Figures/DE_lncRNAset_distribution.pdf")
hist(sig_pairs, xlim = range(0,800))
abline( v = lncRNA_ASD_pairs, col = "red", lwd=2)
dev.off()
