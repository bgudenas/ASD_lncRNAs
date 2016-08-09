# Cyto_converter tests ----------------------------------------------------
source("cyto_converter/cyto_converter.R")
###############################################################
##### Testing cyto_converter and associated functions for accuracy
setwd("../cyto_converter/")



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



# Test Cyto_converter against hg18 bands -----------------------------------------
download.file(url = "http://hgdownload.cse.ucsc.edu/goldenpath/hg18/database/cytoBandIdeo.txt.gz", destfile = "cytoBand_HG18.txt.gz" , mode = "wb")
R.utils:::gunzip("cytoBand_HG18.txt.gz", overwrite = TRUE)
cytobandsHG18 = read.table(file = "cytoBand_HG18.txt", sep = "\t")

# Clean up cytobands ------------------------------------------------------
colnames(cytobandsHG18) = c("Chromosome","Start","End","Band","Stain")
cytobandsHG18$Chromosome = as.character(stringr::str_sub(cytobandsHG18$Chromosome, start=4))
### make band entries  written in full with chromosome prefix Ex. 1q43.1
cytobandsHG18$Band = paste0(cytobandsHG18$Chromosome, cytobandsHG18$Band)
cytobandsHG18 = cytobandsHG18[nchar(cytobandsHG18$Chromosome) <= 2, ] ## remove alt haplotypes
cytobandsHG18 = cytobandsHG18[cytobandsHG18$Chromosome != "M",  ]


table(cytobands$Start ==  cytobandsHG18$Start)
table(cytobands$End ==  cytobandsHG18$End)

test2 = cyto_converter(bands = cytobandsHG18$Band, cytobands = cytobands)
## check against UCSC cytoband file and verify output matches exactly
c(table(test2$Start==cytobands$Start), table(test2$End==cytobands$End))

cyto_converter(bands = "16p11.2", cytobands = cytobandsHG18)
cyto_converter(bands = "15q11-q13", cytobands = cytobandsHG18)




# test against SFARI CNVs -------------------------------------------------
SFARI = read.csv("cnv-summary.csv") ## from SFARI w/ added first column containing row.nums
SFARI_cyto = cyto_converter(SFARI$CNV.Locus, cytobands = cytobands)

unmatched = SFARI_cyto[rowSums(is.na(SFARI_cyto)) > 0, ]
# unmatched
# Cytoband Chromosome    Start End
# 2344    3q26.32-q33          3 1.76e+08  NA
# 4896         8p13.1          8       NA  NA
# 5839         9q34.4          9       NA  NA
# 7326       12q25.32         12       NA  NA
# 10560 22q11.2-q22.3         22 1.74e+07  NA
# manually check bands with NAs
    # 1. 3q26.32-q33: 3q33 does not exist in hg38
    # 2. 8p13.1: 8p13 exceeds boundary of 8p12 in hg38 
    # 3. 9q34.4: exceeds boundary
    # 4. 12q25.32: exceeds boundary
    # 5. 22q11.2-q22.3: end point 22q22.3 exceeds boundary
SFARI_cyto = SFARI_cyto[rowSums(is.na(SFARI_cyto)) == 0, ]
# make sure all CNV widths are positive
table((SFARI_cyto$End - SFARI_cyto$Start)>=0)
TRUE 
11364



