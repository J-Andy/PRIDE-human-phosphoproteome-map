## ---------------------------
##
## Script name: pride-mq-phospho-to-uniprot.R
##
## Purpose of script: process MaxQuant phosphoproteomics result from Ochoa et.al. to include in UniProt
##
## Authors: Dr. Andrew Jarnuczak and Dr. David Ochoa
##
## Date Created: 2019-08-30
##
## Copyright (c) Andrew Jarnuczak, 2019
## Email: jarnuczak@ebi.ac.uk
##
## ---------------------------
##
## Notes: 
## The sctript uses data from Ochoa et. al. (https://www.biorxiv.org/content/10.1101/541656v1) 
## The relevant MQ output files can be downloaded from "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/02/PXD012174/txt-001PTM.zip"
## It appears that in the msmsm.txt table MaxQuant does not report infromation for the phosphosites with Localization probability <0.5. Since we extract the metrics from msms.txt table, there is no Maximum Search Engine Score, PXD, Biological Sample and PUBMEDID information for these sites in the final output. Also the number of Spectral counts is 0 in these cases.
##
## ---------------------------

# load packages
packages <- c("RCurl" ,"data.table", "tidyverse", "plyr")
for(i in packages){
  print(i)
  if(!require(i, character.only = TRUE)){
    install.packages(i, repos="https://cran.ma.imperial.ac.uk/")}
  library(i, character.only = T)
}

# to download original MQ output from PRIDE go to "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2019/02/PXD012174/txt-001PTM.zip"

# there are two metadata files (fileInfo and pxdInfo) we need for annotations (with PRIDE project information, biological sample origin and Pubmed IDs) 
fileInfoURL <- "https://raw.githubusercontent.com/J-Andy/PRIDE-human-phosphoproteome-map/master/metadata/PRIDE-human-phospho-datasets-fileInfo.csv"
fileInfo <- read.csv( text = getURL(fileInfoURL), stringsAsFactors = FALSE) %>%
  separate(Experiment, c("biological_sample", "PXD"), sep = "_")

pxdInfoURL <- "https://raw.githubusercontent.com/J-Andy/PRIDE-human-phosphoproteome-map/master/metadata/PRIDE-human-phospho-datasets-pxdInfo.csv"
pxdInfo <- read.csv(textConnection(getURL(pxdInfoURL)), stringsAsFactors = FALSE)
#
# define txt directory
pridedir <- "E:/phospho-analysis/R-project/txt/"

# read in msms.txt MQ output. Note this is a large file and will take a long time to read
msms <- fread(paste("cut -f  1,2,4,6,10,11,35-40,62,66-70,73 ", pridedir, "msms.txt", sep = ""), verbose = TRUE) %>%
  setNames(make.names(names(.))) %>%
  ## filter(Reverse != "+") %>%
  #    filter(Phospho..STY..site.IDs != "") %>%
  left_join(fileInfo, by = c("Raw.file" = "Name"))
# saveRDS(msms, "msms_annotated.rds")
# rm(msms)
# gc()
# msms <- read_rds("msms_annotated.rds")

# read Phospho (STY)Sites.txt MQ output -- this is where the key phosphosite information is captured
sty  <- fread(paste(pridedir, "Phospho (STY)Sites.txt", sep = ""), showProgress = FALSE, stringsAsFactors = FALSE)
sty <- setNames(sty, make.names(names(sty)))


# combine Phospho (STY)Sites.txt and the annotated msms.txt tables
all_msms <- sty %>%
  select(id, MS.MS.IDs, Score) %>%
  rename(  "STY_table_score" = "Score") %>%
  mutate(MS.MS.IDs = str_split(MS.MS.IDs, ";")) %>%
  unnest()   %>%
  distinct()  %>%
  left_join(msms %>%
              select(id, PXD, biological_sample, Score) %>%
              mutate(id = as.character(id)),
            by = c("MS.MS.IDs" = "id")) %>%
  left_join(pxdInfo %>% select(Project_ID, Pubmed),
            by = c("PXD" = "Project_ID")) %>%
  group_by(id)  %>%
  summarise(Spectralcounts = n(), 
            PXDs = paste(unique(PXD), collapse = ";") ,
            Biological_sample = paste(unique(biological_sample), collapse = ";") ,
            Pubmeds = paste(unique(Pubmed), collapse = ";"), 
            Max_search_score_MSMSid_Table = max( Score), 
            STY_table_score = max(STY_table_score))  # %>% mutate(Spectralcounts = ifelse( PXDs == "NA", 0, Spectralcounts))

# add the most interesting scores, the corresponding datasets in PRIDE, PMIDs if available, etc. 
out <- sty %>%
  filter(Reverse != "+") %>%
  filter(Potential.contaminant != "+") %>%
  select(Protein, Position, Amino.acid, id, Localization.prob, PEP, Proteins, Positions.within.proteins, Position.in.peptide, Best.localization.MS.MS.ID) %>%
  left_join(msms %>%
              select(id, Modified.sequence),
            by = c("Best.localization.MS.MS.ID" = "id")) %>%
  select(-Best.localization.MS.MS.ID) %>%
  left_join(all_msms, by = "id") %>% 
  mutate(Modified.sequence = gsub("\\(..\\)", "", Modified.sequence) )%>% 
  mutate(Modified.sequence = gsub("_", "", Modified.sequence) )

head(out)
colnames(out)

out.aggregated.by.peptide <-  ddply(out, .(Modified.sequence), summarize, 
                                    Protein=paste(Protein, collapse=";"), 
                                    Position=paste(Position, collapse=";"), 
                                    Amino.acid=paste(Amino.acid, collapse=";"), 
                                    id=paste(id, collapse=";"),
                                    Localization.prob=paste(Localization.prob, collapse=";"), 
                                    PEP=paste(PEP, collapse=";"),
                                    Proteins=paste(Proteins, collapse=";"), 
                                    Positions.within.proteins=paste(Positions.within.proteins, collapse=";"),
                                    Position.in.peptide=paste(Position.in.peptide, collapse=";"),
                                    Spectralcounts=paste(Spectralcounts, collapse=";"),
                                    PXDs=paste(PXDs, collapse=";"),
                                    Biological_sample=paste(Biological_sample, collapse=";"),
                                    Pubmeds=paste(Pubmeds, collapse=";"),
                                    Max_search_score_MSMSid_Table=paste(Max_search_score_MSMSid_Table, collapse=";"),
                                    STY_table_score=paste(STY_table_score, collapse=";")
                                    
                                    )


#
write.csv(out, file = "pride-mq-phospho-to-uniprot-out.csv", row.names = FALSE)

write.csv(out.aggregated.by.peptide, file = "pride-mq-phospho-to-uniprot-out-aggregated-by-peptide-seq.csv", row.names = FALSE)
