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

lnc_mods =c("midnightblue","blue","brown","green")
setwd("./Modules/")
for (mod in lnc_mods){
    setwd(mod)
    GO = read.table(paste0(mod,"_raw.txt"), sep = "\t", header = TRUE)
    GO = GO[grepl("GOTERM_BP_FAT",GO$Category) | grepl("GOTERM_CC_FAT",GO$Category) | grepl("GOTERM_MF_FAT",GO$Category), ]
    GO=GO[,-6]
    GO_keep = GO[1:10,c(1,2,12)]
    GO_keep$FDR = -log(GO_keep$FDR)
    
    GO_keep$Term = as.character(GO_keep$Term)
    for ( i in 1:nrow(GO_keep)){
        term =tail(strsplit(as.character(GO_keep$Term[i]), "~")[[1]],1)
        GO_keep$Term[i]=term
    }
    
    GO_keep$color = "blue"
    GO_keep$color[grepl("GOTERM_CC_FAT",GO_keep$Category)] = "red"
    GO_keep$color[grepl("GOTERM_MF_FAT",GO_keep$Category)] = "green"
    
    pdf(paste0(mod,"GO_terms.pdf"))
    par(mar = c(5, 14, 2, 2));
    barplot(GO_keep$FDR, names.arg = GO_keep$Term, las=1, col=GO_keep$color,horiz = TRUE, main=paste(mod, "gene ontology enrichment"), xlab = "-Log(FDR-adjusted P-value)")
    abline(v=3, lwd=3)
    dev.off()
    setwd("..")
    
    }
