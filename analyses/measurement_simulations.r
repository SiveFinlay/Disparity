#04/07/2014

#Script to test my simulations code on linear measurements
  #Data file is my mandible measurements
  #(created by taking the intact mandible measurements from Skulls_after_remeasuring_06_2013
  #and combining it with the mandible measurements from Skulls_FMNH_Sept2013)

library(ape)
library(geiger)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PValueFunction_FromDistribution.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r" )

#---------------------------------------------------------
#phylogenies
setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
    mytrees <- read.tree("Mands_tenrec+gmole_101trees.phy")

#data
setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/mands")

  #measurements
    mands <- read.csv("Mandibles_measurements.csv")

  #taxonomic information
    taxa <- read.csv("Mands_14_03_2014_Images+Specimens.csv", header=T)
    
#------------------------------------------------
#Taxonomic information for data that I have in the measurements

spec.mands <- levels(mands$SpecID)
spec.taxa <- levels(taxa$SpecID)

#find the common species in the two data sets
  #spec ID in spec.taxa but not spec.mands
    rem.taxa <- which(match(spec.taxa,spec.mands, nomatch=1000)==1000)
  #SpecID that needs to be removed from taxa
    rem.taxa.specid <- spec.taxa[rem.taxa]
  #Remove that SpecID from taxa
    taxa.updated <- taxa[-(which(taxa$SpecID==rem.taxa.specid)),]
  #Drop the unused levels
    taxa.updated <- droplevels(taxa.updated)

  
  #spec ID in spec.mands but not spec.taxa
    rem.mands <- which(match(spec.mands,spec.taxa, nomatch=1000)==1000)
  #SpecIDs that need to be removed from mands
    rem.mands.specid <- spec.mands[rem.mands]
  #Find the corresponding row numbers 
    rem.id <- list(NA)
      
      for (i in 1:length(rem.mands.specid)){
        rem.id[[i]] <- which(mands$SpecID==rem.mands.specid[i])
      }
  
  #Remove those specimens from mands
    mands.updated <- mands[-unlist(rem.id),] 
  #Drop the unused levels
    mands.updated <- droplevels(mands.updated)
    
#Check that the two data frames now have the same variables

  spec.mands.updated <- levels(mands.updated$SpecID)
  spec.taxa.updated <- levels(taxa.updated$SpecID)
  
  spec.mands.updated==spec.taxa.updated #all true
  
#-----------------------------------------------------------------
#Combine the two data sets together; measurements from mands.updated and taxonomy from taxa.updated
  #use matching SpecID values to combine

  combine <- merge(mands.updated, taxa.updated)
  
#Select the golden mole and tenrec species only
  tc.gm <- combine[c(which(combine$Family_05 == "Chrysochloridae"), which(combine$Family_05 == "Tenrecidae")),]
  tc.gm <- droplevels(tc.gm)

#Get the median values for each measurement

#List of the different SpecIDs
  specid <- levels(tc.gm$SpecID)

#List of the different measurements
  measures <- levels(tc.gm$Measure)

#Matrix for the median values
  median.list <- matrix(NA,(length(measures)*length(specid)),3)
  colnames(median.list) <- c("SpecID", "Measurement", "Median")

for (i in 1:length(specid)){
  for (j in 1:length(measures)){
    median.list[(j+(i*(length(measures)) - length(measures))),1] <- specid[i]
    median.list[(j+(i*(length(measures)) - length(measures))),2] <- measures[j]
    median.list[(j+(i*(length(measures)) - length(measures))),3] <- median(tc.gm$Value[which(tc.gm$SpecID == specid[i] & test$Measure == measures[j])])    
  }
}

#Get rid of the speech marks around the character objects
  median.list <- as.data.frame(median.list)

#Add the taxonomic information
  median.taxa <- merge(median.list, taxa.updated)
  
#Select only the relevant columns
  mydata <- median.taxa[,c(1,2,3,7,8,9,11,13)]
  
#-------------------------------------------------------------------------------
#Check that the data and the phylogenies have the same species
TreeOnly <- tree.only(mytrees,mydata$Binomial_05)        #null

#species which are in the data but not in the trees
  DataOnly<-as.list(rep(NA,length(mytrees)))

  for (i in 1: length(mytrees)){
    DataOnly[[i]]<-setdiff(mydata$Binomial_05, mytrees[[i]]$tip.label)
  }
  
#Remove these two _sp. species from the data
  rem.sp <- c(which(mydata$Binomial_05 == "Chrysochloris_sp."), which(mydata$Binomial_05 == "Oryzorictes_sp."))
  mydata <- droplevels(mydata[-(rem.sp),])
   
#--------------------------------------------------------
#Simulate character evolution

#One trait on its own: mandible length
  mand.length <- mydata[which(mydata$Measurement=="15_ML"),]
  mand.length <- droplevels(mand.length)

#Unique species names 
 binom<-levels(mand.length$Binomial_05)

#Average values for each species
  mand.length.mean <- matrix(NA, length(binom),1)
  for (i in 1:length(binom)){
    mand.length.mean[i,1] <- mean(as.numeric(as.vector(mand.length$Median[which(mand.length$Binomial_05==binom[i])])))
  }

  #add the species as rownames
  rownames(mand.length.mean) <- binom

#Separate variance covariance matrix of mandible length for each of the phylogenies

#Separate variance covariance matrix of the shape data for each of the phylogenies
  varcov <- as.list(rep(NA,length(mytrees)))
    for(i in 1:length(mytrees)){
      varcov[[i]] <- vcv(phy=mytrees[[i]],mand.length.mean)
    } 
    
#simulate shape evolution on each phylogeny
  length.sim <- as.list(rep(NA,length(mytrees)))
    for (i in 1: length(mytrees)){
      length.sim[[i]] <- sim.char(mytrees[[i]], varcov[[i]],nsim=1000,model="BM")
    }
    
#Combine simulations into one list
  simlist <- list.arrays.to.matrices(length.sim)

#These simulations work but it would make more sense to have more than one trait
  #then I can compare them with PC analyses or some other metric   