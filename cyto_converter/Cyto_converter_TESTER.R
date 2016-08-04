# Cyto_converter tests ----------------------------------------------------
source("cyto_converter/cyto_converter.R")
###############################################################
##### Testing cyto_converter and associated functions for accuracy



## test CNV_ter over all possible combinations of chromosomes qters and pters 
results = c()
 for (i in unique(cytobands$Chromosome)){
     band=paste0(i,"pter")
     #print(paste0(i,"pter"))
     pter = CNV_ter(band, chrom=i)
     
      band=paste0(i,"qter")
      #print(paste0(i,"qter"))
      qter = CNV_ter(band, chrom=i)
     ## test if--------------------
         #all pters start at 0 
         #the end of pters match the start of qters
         #the end of qters match the end of same chromosome
     results = c(results, (pter[1]==0 & pter[2] == qter[1] & (max(cytobands$End[cytobands$Chromosome==i]) == qter[2]) ) ) 
 }
print(table(results))   


# Test Cyto_converter against hg38 bands ----------------------------------
test1 = cyto_converter(bands = cytobands$Band, cytobands = cytobands)
## check against UCSC cytoband file and verify output matches exactly
c(table(test1$Start==cytobands$Start), table(test1$End==cytobands$End))
#  TRUE TRUE



# Test Cyto_converter against hg19 bands -----------------------------------------
download.file(url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/cytoBand.txt.gz", destfile = "cytoBand_HG19.txt.gz" , mode = "wb")
R.utils:::gunzip("cytoBand_HG19.txt.gz")
cytobandsHG19 = read.table(file = "cytoBand_HG19.txt", sep = "\t")

# Clean up cytobands ------------------------------------------------------
colnames(cytobandsHG19) = c("Chromosome","Start","End","Band","Stain")
cytobandsHG19$Chromosome = as.character(stringr::str_sub(cytobands$Chromosome, start=4))
### make band entries  written in full with chromosome prefix Ex. 1q43.1
cytobandsHG19$Band = paste0(cytobands$Chromosome, cytobands$Band)
cytobandsHG19 = cytobands[nchar(cytobands$Chromosome) <= 2, ] ## remove alt haplotypes
cytobandsHG19 = cytobands[cytobands$Chromosome != "M",  ]

test2 = cyto_converter(bands = cytobandsHG19$Band, cytobands = cytobands)
## check against UCSC cytoband file and verify output matches exactly
c(table(test2$Start==cytobands$Start), table(test2$End==cytobands$End))






