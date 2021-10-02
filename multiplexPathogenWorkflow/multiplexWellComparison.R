
# Author: Henry Becker
# Date: 7/13/2021
# Purpose: To streamline the analysis of duplex multiplex DNA pathogen testing results for ticks


#############################
 #          STEPS          #
#############################
# 1. Open Bio-Rad CFX Manager on lab desktop computer
#
# 2. Drag and drop the qPCR file from the flashdrive into the Bio-Rad CFX Manager application
#
# 3. Select the end point data view tab 
#
# 4. Select "custom Export..." from the Export menu at the top of the screen
# / the goal is to export a data table in the correct format for this R script to work with it\
# 
# 5. Near the top of the export menu unclick the box for a header 
# /the R script looks for the column names in the first row of the data file. 
# If the file has a header the script will not be able to determine the start of the columns\
# 
# 6. Make sure the necessary columns are selected for export by clicking the boxes 
# /Necessary columns:
# "Well", "Fluor", "Content", "Sample", "Cq", "Call", "End.RFU"\
# 
# 7. Make sure the file name of the data table does not have any spaces
# /pro tip: file names can never have spaces!\
#
# 8. Export the data table
# /you can export the file to any folder you want but you must set the working directory of 
# this R script to that same folder so the script can find the data file\
# 
# 9. Set the working directory of this script to the same folder to which the data table is located
# Hard way: use the "setwd()" command (you need to know the directory structure of the computer) 
# Easy way: At the top of RStudio select the "Session" tabe -> select "Set Working Directory" -> select "Choose Directory"
# setting the working directory is critical

# 10. Load the date table into R
# Change the name of the green text below to the file name of the data table you exported, then run the two lines below
noHeader <- read.csv(file = "file_name_here.csv")

noHeader
# if you have succesfully loaded the table it print in the console window, 
# if you get the error "Error: object 'noHeader' not found" you've gotta go back through the previous steps
# and make sure you didn't miss anything

# 11. Set the name of the output PDF by changing the run number below (if you don't change this each time you will overwrite the last PDF you created)
pdfName <- "batchXviz"

# 12. Run the rest of the script straight through 

##########################################################
 #          Brief Overview of the workflow below       #
##########################################################
# •load qPCR data table into R 
# •subset data table into each pathogen 
# •define a function which takes the pathogen data table as a parameter 
# •the function creates and populates a 384 well matrix with the relevant positives 
# •define a function that takes each matrix as a parameter
# •the function graphs each matrix and colors the matrix according to the color of the dye in the table
# •finally export graphs as PDFs with the pathogen name in the filename


# load the following packages which are used below

# install them if you havent yet 
#install.packages("stringr")
#install.packages("reshape")
#install.packages("ggplot2")

library(stringr)
library(reshape)
library(ggplot2)

# If you did want to keep the data table header for some reason you'd need to skip the first 19 lines to properly load it into R
# withHeader <- read.csv(file = "qPCRheader.csv", skip = 19)

#split "Well" column into a row (letter) and col (number) column 
noHeader$letters <- (str_extract(noHeader$Well, "[aA-zZ]+"))
noHeader$numbers <- as.numeric((str_extract(noHeader$Well, "[0-9]+")))
noHeader$positive <- ifelse(noHeader$Call == "(+) Positive", TRUE, FALSE)
      
# create a new table with only the borrelia data (borrelia uses the FAM dye)
borrelia <- noHeader[grepl("FAM", noHeader$Fluor), ]

# create a new table with only the anaplasma data (anaplasma uses the HEX dye)
anaplasma <- noHeader[grepl("HEX", noHeader$Fluor), ]

# create a new table with only the babesia data (babesia uses the Cy5 and Texas Red dye)
babesia <- noHeader[grepl("Cy5|Texas Red", noHeader$Fluor), ]

# make a string of letters to represent the rows on the wellplate 
letters <- LETTERS[seq( from = 1, to = 16 )]

########## color conversion function  #########

# this function converts a dye color to a representative number to be placed in a matrix
# inputs: dye color (string)
# outputs: number to represent dye (numeric)

# FAM = 1 = green
# HEX = 2 = orange 
# Texas Red = 3 = red
# Cy5 = 4 = blue

getColor <- function(dye){
  if(dye == "FAM"){
    return(1)
  }
  else if (dye == "HEX"){
    return(2)
  }
  else if (dye == "Texas Red"){
    return(3)
  }
  else if (dye == "Cy5"){
    return(4)
  }
  else{
    print("Wait a minute.. no color found")
  }
}

########## matrix population function  #########

# this function populates a 16 by 24 matrix of zeroes (representation of the well plate) 
# with numbers that correspond to the dye if a cell is positive
# inputs: a data table of a specific pathogen (borrelia, anaplasma, or babesia) (data frame)
# outputs: a 16 by 24 matrix with numbers representing positives (matrix)


matrixPopulator <- function(pathogen) {
  
  # create the matrix
  wellPlate <- matrix(0, nrow = 16, ncol=24)
  
  # add letters as rownames and numbers as colnames
  row.names(wellPlate) <- letters
  colnames(wellPlate) <- 1:24
  
  # iterate through each row of the pathogen data table
  # if a row is positive place a number in the corresponding coordinate of the matrix
  for(i in 1:nrow(pathogen)){
    if(pathogen$positive[[i]]==TRUE){
      row = which(letters == pathogen$letters[[i]])
      col = pathogen$numbers[[i]]
      wellPlate[row, col] = getColor(pathogen$Fluor[[i]])
    }
  }
  return (wellPlate)
}

##################

# run the matrix function on each pathogen data table
BorreliaPlate <- matrixPopulator(borrelia)
AnaplasmaPlate <- matrixPopulator(anaplasma)
BabesiaPlate <- matrixPopulator(babesia)

# view each matrix
BabesiaPlate
AnaplasmaPlate
BorreliaPlate

######################

# set the colors that will be graphed for each dye 
cols <- c("0" = "gray0", "1" = "chartreuse4", "2" = "orange2", "3" = "red4", "4" = "lightseagreen")

###################### graphing function  ###################

# this function graphs a 16 by 24 matrix
# inputs: plate=a 16x24 matrix, bool=true/false which determines wether 
# the cell coordinates are included in the graph or left out

graphPlate <- function(plate, bool){
  
  name = deparse(substitute(plate))
  
  df <- melt.matrix(plate)
  colnames(df) <- c("y", "x", "value")
  
  
  df$y = with(df, factor(y, levels = rev(levels(y))))
  
  
p1 <-  ggplot(df, aes(x = factor(x), y = factor(y), fill = factor(value))) + 
    labs(title = name, caption = "filename here") +
    geom_tile(color = "gray") + 
    coord_fixed() +
    geom_text(aes(label=paste(y,x)), size = 2) +
    scale_fill_manual(values = cols, labels = c("Empty/Negative", "FAM", "HEX", "Texas Red", "Cy5"), name = "Fluorophore") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 14),
          plot.caption = element_text(hjust = 0.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

p2 <- ggplot(df, aes(x = factor(x), y = factor(y), fill = factor(value))) + 
  labs(title = name, caption = "filename here") +
  geom_tile(color = "gray") + 
  coord_fixed() +
  scale_x_discrete(position = 'top') + 
  scale_y_discrete() +
  scale_fill_manual(values = cols, labels = c("Empty/Negative", "FAM", "HEX", "Texas Red", "Cy5"), name = "Fluorophore") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.caption = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())

    if(bool == TRUE ){
        return(p1)
    }
    else if (bool == FALSE){
        return(p2)
    }

}

# FAM = 1 = green
# HEX = 2 = orange 
# Texas Red = 3 = red
# Cy5 = 4 = blue

################################

# turn warnings off because the melt warning in the reshape package does not apply here
#options(warn = -1)

# graph each matrix with and without cell coordinates
graphPlate(BorreliaPlate, TRUE)
graphPlate(BorreliaPlate, FALSE)
graphPlate(AnaplasmaPlate, TRUE)
graphPlate(AnaplasmaPlate, FALSE)
graphPlate(BabesiaPlate, TRUE)
graphPlate(BabesiaPlate, FALSE)

# open a PDF file 

filename <- paste0(pdfName, "Viz.pdf")

pdf(file = filename)

# plot each matrix within the PDF
plot(graphPlate(BorreliaPlate, FALSE))
plot(graphPlate(BorreliaPlate, TRUE))
plot(graphPlate(AnaplasmaPlate, FALSE))
plot(graphPlate(AnaplasmaPlate, TRUE))
plot(graphPlate(BabesiaPlate, FALSE))
plot(graphPlate(BabesiaPlate, TRUE))

# Run dev.off() to close the file!
dev.off()

# set warnings back on 
defaultW <- getOption("warn")
options(warn = defaultW)

# YOUR PDF WITH THE WELL PLATE VISUALIZATIONS WILL BE CREATED IN THE WORKING DIRECTORY FOLDER 
# heres the folder in case you forgot
getwd()









