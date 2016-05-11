# Brian Gudenas
# 5/11/16 
#--- Developmental expression plots

workdir = "C:/Users/Brian/Documents/RNAseq/Autism/WGCNA/"
setwd(workdir)
load(file="./Data/Post_geneinfo_network.RData")

library(ggplot2)

p = qplot(x = clinical$MonthsAge[1:185] , y = MEs0[,3][1:185])
p = p + geom_jitter(alpha=1,shape=1)
p = p + geom_smooth(method=loess,size=1,col="blue",span=1/3, alpha=.5)
p
