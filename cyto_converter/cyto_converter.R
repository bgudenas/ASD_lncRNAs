# Cyto_Converter ----------------------------------------------------------
# Brian Gudenas. 8/1/16
# Purpose: convert cytobands / locus to genomic coordinates

# Libraries ---------------------------------------------------------------
library(R.utils)
library(stringr)

#setwd("C:/Users/Brian/Documents/RNAseq/Autism/cyto_converter/")

# download.file(url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz", destfile = "cytoBand_HG38.txt.gz" , mode = "wb")
# R.utils:::gunzip("cytoBand_HG38.txt.gz")
cytobands = read.table(file = "C:/Users/Brian/Google Drive/ASD_lncRNAs/cyto_converter/cytoBand_HG38.txt", sep = "\t")

# Clean up cytobands ------------------------------------------------------
colnames(cytobands) = c("Chromosome","Start","End","Band","Stain")
cytobands$Chromosome = as.character(stringr::str_sub(cytobands$Chromosome, start=4))
### make band entries are written in full with chromosome prefix Ex. 1q43.1
cytobands$Band = paste0(cytobands$Chromosome, cytobands$Band)
cytobands = cytobands[nchar(cytobands$Chromosome) <= 2, ] ## remove alt haplotypes
cytobands = cytobands[cytobands$Chromosome != "M",  ] ## remove mitochondria


### CNV_ter checks for "qter"s or "pter"s which indicate entire chromosomal arms affected and outputs the entire range of the affected arm
CNV_ter = function(band, chrom){
    if (grepl("pter", band)) {
        e_pos = max(cytobands$End[cytobands$Chromosome == chrom  & grepl(paste0(chrom,"p"),cytobands$Band)])
        s_pos = 0
    }else if (grepl("qter",band)) {
        e_pos = max(cytobands$End[cytobands$Chromosome == chrom & grepl(paste0(chrom, "q"), cytobands$Band)])
        s_pos = min(cytobands$Start[cytobands$Chromosome == chrom & grepl(paste0(chrom, "q"), cytobands$Band)])
}
return(c(s_pos,e_pos))
}

## test CNV_ter over all possible combinations of chromosomes
# for (i in unique(cytobands$chromosome)){
#     band=paste(i,"pter",sep="")
#     print(paste(i,"pter"))
#     print(CNV_ter(band))
#     band=paste(i,"qter",sep="")
#     print(paste(i,"qter"))
#     print(CNV_ter(band))
#     
# }



#### this function deals with bands that are too specific for the reference cytobands EX 1q43.11 becomes 1q43.1
### if a match is not found after degeneration; CNV_ter is tried; else return NA
degenerate = function(band, chrom){
    original = band  ## make a copy of the input
    while (!(band %in% cytobands$Band) & nchar(band)>1) {
        band = substr(band, 1, nchar(band)-1)  ## removes last character from string
        if (substr(band, nchar(band), nchar(band)) == ".") {  ## makes sure it does not end on a period
            band = substr(band,1, nchar(band)-1)            ## or else it removes another character
        }
    }
    if (band %in% cytobands$Band){
        s_pos = cytobands[cytobands$Band == band, 2]
        e_pos = cytobands[cytobands$Band == band, 3]
    }else if (grepl("pter|qter", original)){ ######## try CNV_ter (check for "q|p ters" else leave NAN)
                ter_check= CNV_ter(original, chrom)
                s_pos=ter_check[1]
                e_pos=ter_check[2]
    }else {  
        s_pos= NA
        e_pos = NA
    }
    return( c(s_pos,e_pos))
}

cyto_converter = function(bands, cytobands){
#bands is a list of cytobands for conversion (including chromosome prefix)
#cytobands is a reference dataframe from UCSC with 1st,2nd,3rd,4th columns being Chromosome, Start, End, Band
    bands = as.character(bands)
    starts_vec= vector(mode = "list", length = length(bands))
    ends_vec= vector(mode = "list", length = length(bands))
    chrom_vec= vector(mode = "list", length = length(bands))
    
    for (i in 1:length(bands)){

        chrom= strsplit(bands[i],"p|q")[[1]][1] #split band on p OR q and select first element
        chrom_vec[i]=chrom
        if (grepl("-",bands[i])) {  ##### this chunk deals with joined cytobands. EX. 1q43-1q42
            band_split = strsplit(bands[i],"-")[[1]]
            band_start = band_split[1]
            band_end = paste0(chrom, band_split[2])
            s_pos = cytobands[cytobands$Band == band_start, 2]
            e_pos = cytobands[cytobands$Band == band_end, 3]
            
            if (!(length(s_pos))) {
                if (sum(grepl(band_start, cytobands$Band[cytobands$Chromosome==chrom]))!=0) {
                    matches = cytobands[grepl(band_start, cytobands$Band), ]
                    starts_vec[i] = min(matches[, 2])
                } else  {
                    starts_vec[i]= degenerate(band_start, chrom)[1]
                        }
                                  } else starts_vec[i] = s_pos
            if (!(length(e_pos))) {
                if (sum(grepl(band_end, cytobands$Band[cytobands$Chromosome==chrom]))!=0){
                    matches = cytobands[grepl(band_end, cytobands$Band), ]
                    ends_vec[i] = max(matches[, 3])
                        } else  {
                                ends_vec[i]= degenerate(band_end, chrom)[2]
                                    }
                                } else ends_vec[i] = e_pos
            
        }else if (bands[i] %in% cytobands$Band == FALSE ) {
            matches = cytobands[grepl(bands[i], cytobands$Band), ]
            if (nrow(matches) == 0){
                degen_pos = degenerate(bands[i], chrom) 
                starts_vec[i] = degen_pos[1]
                ends_vec[i] = degen_pos[2]
                             
            } else if (nrow(matches)!=0){
                starts_vec[i] = min(matches[, 2])
                ends_vec[i] = max(matches[, 3])
                    }                                           
        } else {           
            s_pos = cytobands[cytobands$Band == bands[i], 2]
            e_pos = cytobands[cytobands$Band == bands[i], 3]
            if (length(e_pos) & (length(s_pos))) {
                starts_vec[i] = s_pos
                ends_vec[i] = e_pos
                }
            }
    }
results_df <- data.frame(bands, unlist(chrom_vec), unlist(starts_vec), unlist(ends_vec))
colnames(results_df)<-c("Cytoband","Chromosome","Start","End")

return(results_df)
}

top8 = cyto_converter(c("16p11.2",
                 "15q11-13",
                 "15q13.3",
                 "1q21.1",
                 "22q11",
                 "7q11.23",
                 "17q12",
                 "3q29") , cytobands = cytobands)
