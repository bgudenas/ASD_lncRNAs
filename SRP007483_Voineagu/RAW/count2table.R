
### change workdir to where raw count data is
setwd("C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/RAW/Counts")

count=c()
samples=c()
for (file in list.files()){
    sample=strsplit(strsplit(file,"_")[[1]][3],".txt")[[1]]
    samples=c(samples,sample)
    counts=read.table(file)
    genes=counts$V1
    count=c(count,counts$V2)
}
c_table=matrix(data=count,nrow=length(genes),ncol=length(list.files()))
rownames(c_table)=genes
colnames(c_table)=samples

write.csv(c_table,"C:/Users/Brian/Documents/RNAseq/Autism/SRP007483_Voineagu/RAW/Counts_table.csv")
