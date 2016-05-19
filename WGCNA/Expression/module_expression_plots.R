# Brian Gudenas
# 5/11/16 
#--- Developmental expression plots

workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)
load(file="./Data/Post_geneinfo_network.RData")

library(ggplot2)

setwd("./Figures/")
dir.create("Expression Trajectory")
setwd("Expression Trajectory/")

for (mod in 1:ncol(MEs0)){
    modName = names(table(moduleColors))[mod]
    p = qplot(x = clinical$MonthsAge[clinical$MonthsAge <= 22] , y = MEs0[,mod][clinical$MonthsAge <= 22], main = paste(paste0(toupper(substr(modName, 1, 1)), tolower(substring(modName, 2))), " Module Developmental Expression Trajectory"),
        xlab="Post-Conception Developmental Months", ylab = "Module Eigengene (PC1)")
    p = p + geom_jitter(alpha=.5,shape=1) + theme_bw()
    p = p + geom_smooth(method=loess,size=1,col= modName,span=1/3, alpha=.5)
    pdf(paste0(modName, ".pdf"))
    print(p)
    dev.off()
}
