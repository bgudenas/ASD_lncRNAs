# Brian Gudenas
# 5/11/16 
#--- Developmental expression plots

library(ggplot2)
library(gridExtra)

load(file="./WGCNA/Data/Post_geneinfo_network.RData")

setwd("./WGCNA/Figures/")
dir.create("Expression Trajectory")
setwd("Expression Trajectory/")

for (mod in 1:ncol(MEs0)){
    modName = names(table(moduleColors))[mod]
    
    mod_terms = terms_full[terms_full$Module==modName, ]
    mod_terms = mod_terms[order(mod_terms$Pvalue), ]
    term_names = mod_terms$Term[order(mod_terms$Pvalue)]
    
    p1 = ggplot(mod_terms, aes( y=-log10(Pvalue), x=reorder(Term, Pvalue), fill= type)) +
    geom_bar( stat = "identity", position = position_dodge()) + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+
        coord_flip()+
        ylab(expression(-Log[10]~Pvalue)) + xlab( "GO Term") +ggtitle(paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), " Module Functional Enrichment") )+ 
        scale_fill_manual(values=c("#56B4E9",  "#E69F00"))
    
    df = as.data.frame(cbind(Months=clinical$Months.Post.Conception[clinical$Months.Post.Conception <= 22], Expression = MEs0[,mod][clinical$Months.Post.Conception <= 22]) )
    p2= ggplot(df, aes(x = Months, y=Expression), main = paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), " Module Developmental Expression Profile"),
        xlab="Post-Conception Developmental Months", ylab = "Module Eigengene (PC1)") +
            geom_point(shape=1, col ="grey") + theme_bw() +geom_jitter(alpha=0.7, shape=1)+
            geom_smooth(method=loess, size=1, col= modName, alpha=.5)+
        ggtitle(label = paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), " Module Developmental Expression Trajectory") )+
        geom_vline(xintercept =  10, col="red", lwd=1)
    
    grid.arrange(p1, p2, ncol=2)
    pdf(paste0(modName, ".pdf"))
    print(p)
    
    dev.off()
}
