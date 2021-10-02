# Author: Henry Becker
# Date: 7/5/2021
# Purpose: Simple bar charts for SSRP poster session

setwd("/Users/hensanity/Desktop/MMCRI/posterSession/")

library(viridis)
library(ggplot2)
library(dplyr)
library(tidyverse)

# RedCap BMA data
newData <- read.csv("BMAdata.csv")

# Check data integrity and clean  

newData$bburgdorferi[is.na(newData$bburgdorferi)] <- 0

ulst <- lapply(newData, unique)

print(ulst)

str(newData)

# To Do
# replace NA and redo graphs 
# graphs heidi's results for comparison

### Infection Rates Bar Graph ### 

infections <- newData[c(4:7)]

infections$rowsums <- rowSums(infections)

cleanInf <- infections

table(infections$rowsums)

#cleanInf <- infections[complete.cases(infections), ]

# out of 31 infected nymphs some of those co infected, only 6.5% (2) acquired their infections from deer and 0% from mice   

sums <- cleanInf %>% summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)))

x <- (110 - rowSums(sums))
sums$none <- x
sums

cleanInf$Bor.Bab.Ana <- rowSums(cleanInf[,c("bburgdorferi", "bmicroti", "aphagocytophylum")] == "1")

cleanInf$Bor.Bab <- rowSums(cleanInf[,c("bburgdorferi", "bmicroti")] == "1")
      
cleanInf$Bor.Ana <- rowSums(cleanInf[,c("bburgdorferi", "aphagocytophylum")] == "1")
      
cleanInf$Ana.Bab <- rowSums(cleanInf[,c("bmicroti", "aphagocytophylum")] == "1")
      
cleanInf$Pow.Bor.Bab <- rowSums(cleanInf[,c("bburgdorferi", "bmicroti", "pow")] == "1")

# Pow + Borrelia + Babesia
# Borrelia babesia anaplasma 

# Borrelia + babesia 

# Borellia + anaplasma 
# Babesia + anaplasma 

iRates <- data.frame(pathogen =c ("Powassan virus", "Anaplasma.p", "Babesia.m",  "Borelia.b",  "No pathogen"), num = c(1, 10, 11, 27, 61))
iRates$pct <- iRates$num / sum(iRates$num)
iRates$pathogen <- as.factor(iRates$pathogen)
iRates$pct <- lapply(iRates$pct, round, 5)
iRates$pct <- as.numeric(iRates$pct)


ggplot(iRates, aes(x = reorder(pathogen, num), y = num, fill = pathogen, label = scales::percent(pct))) + 
      geom_col(position = 'dodge') + 
      geom_text(position = position_dodge(width = .9), vjust = -0.5, size = 4) + 
      scale_fill_viridis_d() +
      ggtitle("Infection Rates", subtitle = "I. scapularis nymphs collected at WNERR") + 
      ylab("# of nymphs infected") + 
      scale_y_continuous(breaks=seq(0,60,5)) + 
      theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 20, hjust = .5, vjust = .5, face = "plain"), axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 15)) 


### Blood Meal Analysis Bar Graph ###

bma <- data.frame(animal =c("Unidentified", "Deer", "Mouse"), num = c(100, 7, 3))
bma$pct <- bma$num / sum(bma$num)
bma$animal <- as.factor(bma$animal)

bma$pct <- lapply(bma$pct, round, 5)
bma$pct <- as.numeric(bma$pct)

ggplot(bma, aes(x = reorder(animal, -num), y = num, fill = animal, label = scales::percent(pct))) + 
      geom_col(position = 'dodge') + 
      geom_text(position = position_dodge(width = .9), vjust = -0.5, size = 4) + 
      scale_fill_viridis_d(option = "C") +
      scale_fill_manual(values = c("navy", "yellow", "pink")) +
      scale_fill_manual(values = c("yellow", "pink", "navy")) +
      #scale_fill_manual(values = c("pink", "navy", "yellow")) +
      ggtitle("Blood Meal Analysis", subtitle = "I. scapularis nymphs collected at WNERR") + 
      ylab("# of positive identifications") + 
      scale_y_continuous(breaks=seq(0,100,5)) + 
      theme(axis.text.x = element_text(color = "grey20", size = 13, angle = 20, hjust = .5, vjust = .5, face = "plain"), axis.title.x = element_blank(), legend.position = "none", plot.title = element_text(hjust = 0.5, size = 20), plot.subtitle = element_text(hjust = 0.5, size = 15)) 



### STACKED BAR PLOT ###
#### Deer and mouse contribution to infection rates ###

contribution <- data.frame(blood = c("Unidentified", "Deer"), 
                       Proportion = c(93.5, 6.5), num = c(29, 2))

contribution <- contribution %>% 
  mutate(X = " ")

ggplot(contribution, aes(x = X, y = Proportion, fill = blood)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(num, " Nymphs")),
            position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size = 16) +
  ylab("Percentage") +
  scale_y_continuous(breaks=seq(0,100,10)) +
  xlab(" ") + 
  ggtitle("Proportion of Infections ", subtitle = "By Blood Source") +
  #scale_fill_viridis_d(option = "F", name = "Blood Source", alpha = 1) +
  scale_fill_discrete(name = "Blood Source") +
  scale_x_discrete(breaks = NULL) +
  scale_fill_manual(values = c("red", "orange"), "Blood Source") +
  theme( plot.title = element_text(hjust = .5, size = 18), plot.subtitle = element_text(hjust = .5, size = 15))


### PIE CHARTS ###

# 3D Exploded Pie Chart
library(plotrix)
slices <- c(11, 10, 27, 1, 60) 
lbls <- c("B.microti", "A.phagocytophilum", "B.burgdorferi", "Powassan Virus", "No Pathogen")
pie3D(slices,labels=lbls,
      main="Infection Percentages", theta = 2)

# Pie Chart with Percentages
slices <- c(1, 10, 11, 27, 60) 
lbls <- c("Powassan Virus", "A.phagocytophilum", "B.microti", "B.burgdorferi", "Unidentified")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="I. scapularis Infection Prevalence")

slices <- c(3, 7, 100) 
lbls <- c("White Footed Mouse", "White Tailed Deer", "No Bloodmeal Detected")
pct <- round(slices/sum(slices)*100)
lbls <- paste(lbls, pct) # add percents to labels 
lbls <- paste(lbls,"%",sep="") # ad % to labels 
pie(slices,labels = lbls, col=rainbow(length(lbls)),
    main="Bloodmeal Remnant Identification")

# pie chart of the percent of each infection 

