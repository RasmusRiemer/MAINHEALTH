library(phyloseq)
library(ggplot2)      # graphics
library(readxl)       # necessary to import the data from Excel file
library(dplyr)        # filter and reformat data frames
library(stringr)
library(vegan)
library(ggpubr)
library(metagenomeSeq)
library(tidyr)
library(RColorBrewer)
library(reshape2)
library(utils)
library(metagMisc)


#Set working directory to r-script directory location - will not be run when script is run through source command

#setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
options(getClass.msg=FALSE)

######################################Construct Phyloseq objects##########################################

########### Load raw data

otu_mat <- read.delim2('Data/count_matrix_merged.tsv', skip =0, sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE)
map_mat <- read.delim2("Data/merged_metadata.txt", sep = '\t', stringsAsFactors = FALSE, check.names = FALSE, row.names = 1, header = TRUE, fileEncoding="Latin1") %>%
  type.convert(as.is = TRUE)
str(map_mat)

tax_mat <- read.delim2('Data/taxonomy.tsv', skip =0,
                       sep = '\t', stringsAsFactors = FALSE,
                       check.names = FALSE, row.names = 1,
                       header = TRUE
                       )

# tax_mat <- read.delim2('Data/taxonomy.tsv', skip =0,
#                        sep = '\t', stringsAsFactors = FALSE,
#                        check.names = FALSE, row.names = 1,
#                        header = FALSE
# )
# 
# 
# rownames <- rownames(tax_mat)
# tax_mat[8] <- paste(tax_mat[,1],tax_mat[,2],tax_mat[,3],tax_mat[,4],tax_mat[,5],tax_mat[,6],tax_mat[,7],sep=";") 
# tax_mat <- tax_mat[8]
# head(tax_mat)

tree <- read_tree("Data/tree.nwk")

############### Set 

#str(map_mat)
map_mat$Run <- as.factor(map_mat$Run)
map_mat$Set <- as.factor(map_mat$Set)

############## Format taxonomy

#Remove "k__", "p__" etc from taxonomy column

tax_mat[,1] <- str_remove_all(tax_mat[,1], str_c(c("k__", "p__", "c__", "o__", "f__", "g__", "s__","\\[","\\]"), collapse="|"))


#Split taxonomy column by ";"

split <- colsplit(tax_mat[,1], ";", c("Kingdom","Phylum","Class","Order","Family","Genus","Species"))

tax_mat <- cbind(tax_mat,split)

#Remove old taxonomy column
tax_mat <- tax_mat[,-c(1)]

tax_mat <- as.matrix.data.frame(tax_mat)

otu_mat[] <- mutate_all(otu_mat, function(x) as.numeric(as.character(x)))
#remove empty samples
otu_mat <- otu_mat %>% select_if(colSums(.) != 0)

otu_mat <- as.matrix(otu_mat)

##Construct individual tables is phyloseq format

OTU <- phyloseq::otu_table(otu_mat, taxa_are_rows = TRUE)
TAX <- phyloseq::tax_table(tax_mat)
samples <- phyloseq::sample_data(map_mat)

##Combine in phyloseq object

PSB.all <- phyloseq(OTU, TAX, samples, tree)

PSB.CSS.all <- phyloseq_transform_css(PSB.all)

#Remove left over objects

rm(tax_mat, map_mat, otu_mat, samples, split, OTU, TAX, tree)


#Check mapping file

