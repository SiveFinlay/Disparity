#Test disparity code with anolis data
#Data from dryad repositories of Harmon et al 2010 and Thomas et al 2009
setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data")

source("C:Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")

anole.harmon <- read.table("anolisraw.dat", header=T, sep = "\t", na.strings=".")


anole.thomas <- read.table("Thomas2009_anolisdata.txt", header=TRUE, sep="\t")
anole.thomas.species <- split.binom(as.character(anole.thomas$Species))

#list of species in each data
an.th <- unique(anole.thomas.species$Species)

an.har <- unique(as.character(anole.harmon$NAME))

#find the common species in the two data sets
common <- common.character(an.th, an.har)

