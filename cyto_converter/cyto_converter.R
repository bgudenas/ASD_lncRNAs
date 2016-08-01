rm(list=ls())
setwd("ASD")
rawCNVs <-read.csv("Raw/SFARI/individual-data.csv") ##7/24
cytobands=read.csv("C:/Users//Brian/Documents/ASD/hg38_cytoBand.csv") ## http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoband.txt.gz
##cytoband input is a csv file which has "chromosome","start","end","band" column headers
### band entries are written in full with chromosome prefix Ex. 1q43.1

### CNV_ter checks for "qter"s or "pter"s which indicate entire chromosomal arms affected and outputs the entire range of the affected arm
CNV_ter = function(band){
    chromosome = strsplit(band,"p|q")[[1]][1]
    if (grepl("pter",band)){
        e_pos=max(cytobands$end[cytobands$chromosome == chromosome  & grepl(paste(chromosome,"p",sep=""),cytobands$band)])
        s_pos=0
        return(c(s_pos,e_pos))
    }
    else if (grepl("qter",band)){
        e_pos=max(cytobands$end[cytobands$chromosome == chromosome & grepl(paste(chromosome,"q",sep=""),cytobands$band)])
        s_pos=min(cytobands$start[cytobands$chromosome == chromosome & grepl(paste(chromosome,"q",sep=""),cytobands$band)])
        return(c(s_pos,e_pos))
}
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
### if a match is not found after degeneration; CNV_ter is tried; else return "NAN"
degenerate = function(band){
    original=band  ## make a copy of the input
    while (!(band %in% cytobands$band) & nchar(band)>1){
        band =substr(band,1,nchar(band)-1)  ## removes last character from string
        if (substr(band,nchar(band),nchar(band))==".") {  ## makes sure it does not end on a period
            band =substr(band,1,nchar(band)-1)            ## or else it removes another character
        }
    }
    if (band %in% cytobands$band){
        s_pos=cytobands[cytobands$band==band,2]
        e_pos=cytobands[cytobands$band==band,3]
    }else if (grepl("pter|qter",original)){ ######## try CNV_ter (check for "q|p ters" else leave NAN)
                ter_check= CNV_ter(original)
                s_pos=ter_check[1]
                e_pos=ter_check[2]
    }else {  
        s_pos=as.numeric("NAN")
        e_pos = as.numeric("NAN")
    }
    return(list(s_pos,e_pos))
}

cyto_converter = function(cyto_bands){  ##input is a list of cytobands
    starts_vec= vector(mode = "list", length = length(cyto_bands))
    ends_vec= vector(mode = "list", length = length(cyto_bands))
    chrom_vec= vector(mode = "list", length = length(cyto_bands))
    cyto_bands=as.character(cyto_bands)
    for (i in 1:length(cyto_bands)){
        s_pos="NaN"
        e_pos="NaN"
        chrom=strsplit(cyto_bands[i],"p|q")[[1]][1]
        chrom_vec[i]=chrom
        if (grepl("-",cyto_bands[i])){  ############# this chunk deals with joined cytobands. EX. 1q43-1q42
            start=strsplit(cyto_bands[i],"-")[[1]][1]
            chrom=strsplit(cyto_bands[i],"p|q")[[1]][1]
            end=paste(chrom,strsplit(cyto_bands[i],"-")[[1]][2],sep="")
            s_pos=cytobands[cytobands$band==start,2]
            e_pos=cytobands[cytobands$band==end,3]
            
            if (!(length(s_pos))) {
                if (sum(grepl(start,cytobands$band[cytobands$chromosome==chrom]))!=0){
                    matches= cytobands[grepl(start,cytobands$band),]
                    starts_vec[i] = as.numeric(as.character(min(matches[,2])))
                } else  {
                    degen_pos=degenerate(start)
                    starts_vec[i]=as.numeric(degen_pos[1])
                        }
                                  } else {starts_vec[i]=as.numeric(as.character(s_pos))}
                          
            if (!(length(e_pos))) {
                if (sum(grepl(end,cytobands$band[cytobands$chromosome==chrom]))!=0){
                    matches=cytobands[grepl(end,cytobands$band),]
                    ends_vec[i]=as.numeric(as.character(max(matches[,3])))
                        } else      {
                    
                                degen_pos=degenerate(end)
                                ends_vec[i]=as.numeric(degen_pos[2])
                                    }
                                } else {ends_vec[i] =(as.numeric(as.character(e_pos)))}                               
        }else if (cyto_bands[i] %in% cytobands$band == FALSE & !(grepl("-",cyto_bands[i])) ){
            matches = cytobands[grepl(cyto_bands[i],cytobands$band[cytobands$chromosome==chrom]),]
            if (nrow(matches)== 0){
                degen_pos = degenerate(cyto_bands[i]) 
                starts_vec[i] = as.numeric(degen_pos[1])
                ends_vec[i] = as.numeric(degen_pos[2])
                            } 
            if (nrow(matches)!=0){
                ends_vec[i]=as.numeric(as.character(max(matches[,3])))
                starts_vec[i] = as.numeric(as.character(min(matches[,2])))
                    }                                           
        } else {           
            s_pos=cytobands[cytobands$band==cyto_bands[i],2]
            e_pos=cytobands[cytobands$band==cyto_bands[i],3]
            if (length(e_pos) & (length(s_pos))){
                starts_vec[i]=as.numeric(as.character(s_pos))
                ends_vec[i] = as.numeric(as.character(e_pos))##end
                }
            }
    }
results_df <- data.frame(cyto_bands,unlist(chrom_vec),unlist(starts_vec),unlist(ends_vec))
colnames(results_df)<-c("cytoband","chromosome","start","end")

# WHAT IS THE CODE BELOW DOING?? ------------------------------------------
for (i in nrow(results_df)){
    if (is.numeric(results_df$start[i])==TRUE & is.numeric(results_df$end[i]==TRUE)){
    if (results_df$start[i] > results_df$end[i]){
        results_df[i,2:4]=NA
    }
    }
}
return(results_df)
}
###############################################################
##### this code below is for testing
test1 = cyto_converter(cytobands$band)
## check against UCSC cytoband file and verify output matches exactly
# table(test1$start==cytobands$start)
# table(test1$end==cytobands$end)
#  TRUE TRUE
######################################################\
# form testing script for cyto_converter (all Pters start 0; all qters end at chrom end etc ...)

SFARI=cyto_converter(rawCNVs$CNV.Locus)

##now extract useful CNV metadata columns for future analysis
# colnames(rawCNVs)[c(1,2,5,6,7,8,10,14,19)]
# [1] "chromosome"        "CNV.Locus"         "Case.Control"      "Patient.Age"       "Patient.Gender"    "Primary.Diagnosis"
# [7] "Cognitive.Profile" "CNV.Type"          "Inheritence"      

CNV_meta=rawCNVs[,c(1,2,5,6,7,8,10,14,19)]
CNV_meta=CNV_meta[!rowSums(is.na(SFARI))>0,]
SFARI=SFARI[!rowSums(is.na(SFARI))>0,]   ### remove any rows with atleast 1 NAN value
#  table(SFARI$cytoband==CNV_meta$CNV.Locus)
# TRUE 
# 56845
write.csv(SFARI,file="CNV/SFARI_converted_CNVs_hg38.csv",row.names=FALSE)


