#09/05/14

#Generic script to create phylogenies for the disparity analysis
  #Can be used with any of the data sources (skdors, sklat, skvent, mandibles)
  #Uses files created by the twofamily_disparity script -> in the output folder

#Steps
#1) Read in data
#2) Prune trees partially (remove genera that aren't in the data)
#3) Add missing species to the trees
#4) Remove species that are only in the tree
#5) Create a matrix of phylogenetic distances among species
#6) Create output files
      #101 phylogenies (mytrees0
      #101 phylogenetic distance matrics (phydists)
      
#NB: Change data input and output files depending on the data used

#---------------------------------------------------------------
library(ape)
library(geiger)
library(phytools)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("http://www.phytools.org/add.species.to.genus/v0.1/add.species.to.genus.R")

#distribution of mammal supertrees with resolved polytomies from Kuhn et al 2011
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data")
  supertree<-read.nexus("FritzTree.rs200k.100trees.tre")

############################################
#READ IN THE DATA; different files for each analysis
##########################################
#Skdors
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skdors")

#1) Mean shape coordinates for each species
#  sps.mean <- dget(file="SkDors_tenrec+gmole_sps.mean.txt")
#2) Taxonomy for those coordinates
#  sps.tax <- read.table(file="SkDors_tenrec+gmole_sps.mean_taxonomy.txt")

#-------------------------------------------------

#Sklat
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/sklat")

#1) Mean shape coordinates for each species
#  sps.mean <- dget(file="SkLat_tenrec+gmole_sps.mean.txt")
#2) Taxonomy for those coordinates
#  sps.tax <- read.table(file="SkLat_tenrec+gmole_sps.mean_taxonomy.txt")

#--------------------------------------------------------
#Skvent
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skvent")

#1) Mean shape coordinates for each species
#  sps.mean <- dget(file="SkVent_tenrec+gmole_sps.mean.txt")
#2) Taxonomy for those coordinates
#  sps.tax <- read.table(file="SkVent_tenrec+gmole_sps.mean_taxonomy.txt")

#----------------------------------------
#Mands
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")

#1) Mean shape coordinates for each species
  sps.mean <- dget(file="Mands_tenrec+gmole_sps.mean.txt")
#2) Taxonomy for those coordinates
  sps.tax <- read.table(file="Mands_tenrec+gmole_sps.mean_taxonomy.txt")

###############################################

#Extract all species in the trees
  all.species <- supertree[1]$tip.label[[1]]

#Matrix of genus and species names of all of the tip labels
  sup.spe <- split.binom(all.species)

#Matrix of genus and species names for species in mydata
  morph.spe <- split.binom(as.character(sps.mean$Binom))

###############################################
#PRUNE TREES PARTIALLY
###############################################

#Genera in the morphological data that are currently in the phylogenies
  common <- intersect(sup.spe$Genus, morph.spe$Genus)

#Select just those genera from the supertrees
  com.gen <- as.list(rep(NA, length(common)))
    for (i in 1: length(common)){
      com.gen[[i]] <- which(sup.spe$Genus == common[i])
    }
  
  com.id <- unlist(com.gen)

#list of id numbers that are not common 
  lose.tip <- all.species[-com.id]
  
#Prune the trees down to a more manageable size
  sup.prun <- remove.missing.species.tree(supertree, lose.tip)

###############################################
#ADD MISSING SPECIES TO THE TREES
###############################################

#Microgale in the super trees and in the morphological data
  sup.mic <- sup.spe[which(sup.spe$Genus == "Microgale"),]
  morph.mic <- morph.spe[which(morph.spe$Genus == "Microgale"),]

#Microgale missing from the trees
  miss.mic <- morph.mic[-(which(morph.mic[,2] %in% sup.mic[,2])),]
  microgale <- paste(miss.mic$Genus, miss.mic$Species, sep="_")

#add the missing microgale species at random to the rest of the microgale genus
  sup.prun.micro <- sup.prun
    for(i in 1:length(sup.prun.micro)){
      for (m in 1:length(microgale)){
        sup.prun.micro[[i]] <- add.species.to.genus(sup.prun.micro[[i]],microgale[m],where="random")
      }
    }

#Find the missing golden mole species

#species in the pruned trees
  prun.spe <- split.binom(sup.prun.micro[[1]]$tip.label)

#golden moles in the morphological data
  morph.gm <- subset.matrix(sps.tax, sps.tax$Family, "Chrysochloridae")
  morph.gm.spe <- split.binom(as.character(morph.gm[,2]))


#golden moles missing from the trees
  miss.gm <- morph.gm.spe[-(which(morph.gm.spe[,2] %in% prun.spe[,2])),]
  miss.gm <- paste(miss.gm$Genus, miss.gm$Species, sep="_")

#golden moles currently in the trees
  tree.gm <- setdiff(morph.gm$Binomial, miss.gm)

#Add the first missing golden mole species
  sup.prun.micro.gm <- add.species.to.MRCA(sup.prun.micro, tree.gm, miss.gm[1])

#Add the second missing golden mole species
  sup.prun.micro.gm2 <- add.species.to.MRCA(sup.prun.micro.gm, tree.gm, miss.gm[2])
  class(sup.prun.micro.gm2) <- "multiPhylo"

#Resolve polytomies at random
  sup.prun.micro.gm2.res <- calc.each.list(sup.prun.micro.gm2, multi2di)
  
###############################################
#REMOVE SPECIES THAT ARE ONLY IN THE TREES
###############################################
  TreeOnly <- tree.only(sup.prun.micro.gm2.res, sps.mean$Binom)

#remove species which are only in the trees and not the data
  #remove the same set of missing species from each phylogeny so use TreeOnly[[1]]
  mytrees <- remove.missing.species.tree(sup.prun.micro.gm2.res, TreeOnly[[1]])

#trees and data now have the same species
  #data.only <- setdiff(sps.mean$Binom, mytrees[[1]]$tip.label)

###############################################
#CREATE A MATRIX OF PHYLOGENETIC DISTANCES
###############################################
#Phylogenetic distance matrix for each of the 101 trees
  phydists <- calc.each.list(mytrees, cophenetic)

###############################################
#CREATE OUTPUT FILES
###############################################

#change the directory 
setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")

#Skdors
#1) 101 fully resolved phylogenies with the species from the shape data
#  write.tree(mytrees, file="SkDors_tenrec+gmole_101trees.phy", append=FALSE)
#2) Phylogenetic distance matrices
#  dput(phydists,file="SkDors_tenrec+gmole_101phydists.txt")

#---------------------------------------------
#Sklat
#1) 101 fully resolved phylogenies with the species from the shape data
#  write.tree(mytrees, file="SkLat_tenrec+gmole_101trees.phy", append=FALSE)
#2) Phylogenetic distance matrices
#  dput(phydists,file="SkLat_tenrec+gmole_101phydists.txt")

#---------------------------------------------------
#Skvent
#1) 101 fully resolved phylogenies with the species from the shape data
#  write.tree(mytrees, file="SkVent_tenrec+gmole_101trees.phy", append=FALSE)
#2) Phylogenetic distance matrices
#  dput(phydists,file="SkVent_tenrec+gmole_101phydists.txt")

#--------------------------------------------------------
#Mands
#1) 101 fully resolved phylogenies with the species from the shape data
  write.tree(mytrees, file="Mands_tenrec+gmole_101trees.phy", append=FALSE)
#2) Phylogenetic distance matrices
  dput(phydists,file="Mands_tenrec+gmole_101phydists.txt")