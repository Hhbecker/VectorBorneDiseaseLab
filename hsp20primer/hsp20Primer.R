# Author: Henry Becker
# Date: 6/15/2021
# Purpose: Design new HSP20 primer. 
# Project: Robich Diapause Pilot Study 
# Overview: Positive control failed using HSP20 primer from Busby et al. Pull all HSP20 sequences from I. scapularis 
# records on GenBank using RENTREZ package. Locate the HSP20 region targeted by Busby et al. and widen the target region by ~200 base pairs.

# UNDER CONSTRUCTION

#set workding directory 
setwd("/Users/hensanity/desktop/hsp20Sequences/")

#install and load necessary packages
#install.packages("rentrez")
#install.packages("openPrimeR")
#BiocManager::install("openPrimeR")

library(seqinr)
library(rentrez)
library(tidyr)
library(dplyr)
library(stringr)


# OpenPrimeR documentation
# https://www.bioconductor.org/packages/release/bioc/manuals/openPrimeR/man/openPrimeR.pdf

# Primer design paper
# https://iubmb.onlinelibrary.wiley.com/doi/full/10.1002/bmb.20461


# convert text to fasta 

# 1. search for hsp 20 sequences 
# 2. download mRNA scripts 
# 3. convert mRNA to its DNA equivalent
# 4. convert to complementary DNA 
# 5. search the cDNA for the primers 

# You have primers on the mRNA that reverse transcriptase uses to sequence DNA from the mRNA (cDNA) and then 
# that cDNA somehow becomes double stranded and thats when sybr green binds 

#Primer Length
# The optimal length of primers is generally accepted as 18–24 bp in length. 
# Longer primers will take longer to hybridize, longer to extend, and longer 
# to remove thus producing less amplicons.

### DOWNLOAD MRNA VS DNA

# 1. get links to vectorBase automatically through Rentrez 
# 1.5 make sure you got the right vectorBase link from each NCBI link
# 2. combine individual fasta files of primer sequences into one fasta file 
# 3. align sequences with primers? 
# 4. convert primers to fasta and align them 
# 5. convert DNA 


# first and second primer in paper are forward and reverse

# Rentrez documentation: https://cran.r-project.org/web/packages/rentrez/vignettes/rentrez_tutorial.html

# https://www.ncbi.nlm.nih.gov/nuccore/183936837?&withparts=on&expand-gaps=on
# https://www.ncbi.nlm.nih.gov/nuccore/DS683839 
# https://www.ncbi.nlm.nih.gov/nuccore/DS643745
# https://www.ncbi.nlm.nih.gov/nuccore/DS979611
# https://www.ncbi.nlm.nih.gov/nuccore/DS646824

# downloaded from vectorBase as CDS (Coding Sequences - untranslated regions and introns excluded)


# returns list of all databases searchable through the rentrez package
entrez_dbs()
# returns a summary of the NCBI nucleotide database
entrez_db_summary("nucleotide")
# returns list of all searchable terms for the NCBI nucleotide database
entrez_db_searchable("nucleotide")


hsp20Search <- entrez_search(db = "nucleotide", term = "Ixodes scapularis[ORGN] AND HSP20", retmax=300)

hsp20Summary <- entrez_summary(db="nucleotide", id=hsp20Search$ids)

#extract summaries of every result from above search
titles <- as.data.frame(t(extract_from_esummary(hsp20Summary, c("title", "caption"))))

# wrangle data frame into tidy data 
titles$id <- as.numeric(rownames(titles))
rownames(titles) <- 1:length(titles$title)
titles$title <- as.character(titles$title)
titles$slen <- as.character(titles$slen)
titles$biomol <- as.character(titles$biomol)
#check data type of each column
sapply(titles, class)

#fetch actual sequences in fasta format
all_recs <- entrez_fetch(db="nucleotide", id=titles$id, rettype = "xml", parsed = TRUE)
print(all_recs)

class(all_recs)

forward <- "TAATACGACTCACTATAGGGTACTATGGCACTGTTCCCACTTT"
nchar(sequences = c(forward, reverse), as.string = TRUE, names = c("I. scap HSP20 forward primer", "I. scap HSP20 reverse primer"), file.out = )
f.fasta <- write.fasta()

reverse <- "TAATACGACTCACTATAGGGTACTTCACTTGCTCTTGACAGCA"
nchar(reverse)

# write(all_recs, file="allSequences.fasta")

#make object out of sequence file
file1 <- readDNAStringSet("file1.fasta")

file2 <- readDNAStringSet("file2.fasta")

file3 <- readDNAStringSet("file3.fasta")

file4 <- readDNAStringSet("file4.fasta")

file5 <- readDNAStringSet("file5.fasta")

allSequences <- c(file1, file2, file3, file4, file5)

alignment <- msa(inputSeqs = allSequences, method = "ClustalW", type = "dna", order = "aligned", verbose = TRUE)


msaPrettyPrint(alignment, output="pdf", showNames="none",
               showLogo="none", askForOverwrite=FALSE, verbose=FALSE)


print(allSequences)

