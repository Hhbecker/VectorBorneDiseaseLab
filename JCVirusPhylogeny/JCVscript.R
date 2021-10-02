# Author: Henry Becker
# Date: 7/29/2021
# Purpose: Update Jamestown Canyon Virus Phylogeny based on Libby's new sequence 
# Project: Libby's Jamestown Canyon Virus Paper 
# Overview: Make the most robust and up to date JCV phylogeny by pulling all JCV sequences from NCBI using Rentrez. Align sequences to dereplicate them. 

### BACKGROUND ### 
# Jamestown Canyon Virus is a relatively rare Orthobunyavirus implicated in several human encephalitis cases.
# 
# Family Orthobunyavirus
# Genus Bunyaviridae 
# 
# Genome is comprised of 3 negative sense single strand RNA 
# 
# Small segment codes for nucleocapsid
# Medium segment codes for surface glycoprotein precursor 
# Large segment codes for RNA polymerase  
# 
# Mosquito vector
# 
# Likely several animal reservoirs including migratory birds and deer.

#######

# UNDER CONSTRUCTION


#set workding directory 
setwd("/Users/hensanity/desktop/JCVproject/")

#install and load necessary packages
#install.packages("rentrez")
#install.packages("tidyr")
#install.packages("dplyr")
#install.packages("stringr")
library(rentrez)
library(tidyr)
library(dplyr)
library(stringr)
library(msa)

# Rentrez documentation: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

## Confirm VBDL sequence is small segment

# returns list of all databases searchable through the rentrez package
entrez_dbs()
# returns a summary of the NCBI nucleotide database
entrez_db_summary("nucleotide")
# returns list of all searchable terms for the NCBI nucleotide database
entrez_db_searchable("nucleotide")


JCVsearch <- entrez_search(db = "nucleotide", term = "Jamestown Canyon Virus[ORGN]", retmax=300)

JCVsummary <- entrez_summary(db="nucleotide", id=JCVsearch$ids)

#extract summaries of every result from above search
sequenceTitles <- as.data.frame(t(extract_from_esummary(JCVsummary, c("title", "slen", "biomol"))))

# wrangle data frame into tidy data 
sequenceTitles$id <- as.numeric(rownames(sequenceTitles))
rownames(sequenceTitles) <- 1:length(sequenceTitles$title)
sequenceTitles$title <- as.character(sequenceTitles$title)
sequenceTitles$slen <- as.character(sequenceTitles$slen)
sequenceTitles$biomol <- as.character(sequenceTitles$biomol)
#check data type of each column
sapply(sequenceTitles, class)

# extract rows containing "segment S", "nucleoprotein", or "nucleocapsid"
sSequences <- sequenceTitles[grepl("segment S|nucleoprotein|nucleocapsid", sequenceTitles$title), ]

#fetch actual sequences in fasta format
all_recs <- entrez_fetch(db="nuccore", id=sSequences$id, rettype="fasta")

# write(all_recs, file="allSequences.fasta")

#make object out of sequence file
mySequenceFile <- system.file("exampleAA.fasta", package="msa")

# VBDL JCV sequence has been added to the fasta file (its the first sequence in the file)
mySequences <- readDNAStringSet("allSequences.fasta")

print(mySequences)

firstAlignment <- msa(inputSeqs = mySequences, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)










#align sequences using msa package ClustalW algorithm
#myFirstAlignment <- msa(inputSeqs = "allSequences.fasta", method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)

msaPrettyPrint(myFirstAlignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)

# try making each tree type from the influenza tutorial 
# try aligning them (maybe remove sequence with Ns)

# include geolocation in entrez search 
# possibly get linked publications to find geolocation of each sequence to do a geomap of sequences 



# Distance-based methods
# Just to recap on some background mentioned earlier: Distance-based trees are produced by calculating 
# the genetic distances between pairs of taxa, followed by hierarchical clustering that creates the 
# actual “tree” look. While there are tons of algorithms to choose from when computing distances, 
# there are two popular clustering methods that are used most frequently. Neighbor-joining- taking 
# the two closest nodes of the tree and defines them as neighbors; you keep doing this until all of 
# the nodes have been paired together.
















