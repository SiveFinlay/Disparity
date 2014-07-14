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
    median.list[(j+(i*(length(measures)) - length(measures))),3] <- median(tc.gm$Value[which(tc.gm$SpecID == specid[i] & tc.gm$Measure == measures[j])])    
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

#Unique species names 
 binom<-levels(mydata$Binomial_05)

#Average values for each species
  mydata.mean <- matrix(NA, length(binom)*length(measures),3)
  
  for (i in 1:length(binom)){
    for (j in 1:length(measures)){
      mydata.mean[(j+(i*(length(measures)) - length(measures))),1] <- binom[i]
      mydata.mean[(j+(i*(length(measures)) - length(measures))),2] <- measures[j]
      mydata.mean[(j+(i*(length(measures)) - length(measures))),3] <- mean(as.numeric(as.vector(mydata$Median[which(mydata$Binomial_05 == binom[i] & mydata$Measure == measures[j])])))
    }
  }                                                                   


  #add the species as rownames
  mydata.mean <- as.data.frame(mydata.mean)
  colnames(mydata.mean) <- c("Binom","Measure", "Mean")

#Re-shape the matrix so that each trait is in a separate column
   mydata.mean.rs <- matrix(NA, length(binom), length(measures))
   
   for (i in 1:length(levels(mydata.mean$Measure))){
    rownames(mydata.mean.rs) <- binom
    colnames(mydata.mean.rs) <- levels(mydata.mean$Measure)
    mydata.mean.rs[,i] <- mydata.mean$Mean[which(mydata.mean$Measure == (levels(mydata.mean$Measure)[i]))]
   }

#-------------------------------------------------------
#Variance covariance matrix for the trait data
  #NB: advice from Adam Algar at BES Macro Nottingham: use ic.sigma
      #but that's been deprecated -> use vcv.phylo instead (previously I used vcv)
      #Here using vcv and vcv.phylo give the same results

#One phylogeny first
 one.tree <- mytrees[[1]]
 
 varcov.phylo <- vcv.phylo(one.tree, mydata.mean.rs)
  varcov <- vcv(one.tree, mydata.mean.rs)
#-----------------------------------------------------------
#Separate variance covariance matrix for each of the phylogenies

#Separate variance covariance matrix of the shape data for each of the phylogenies
  varcov.list <- as.list(rep(NA,length(mytrees)))
    
    for(i in 1:length(mytrees)){
      varcov.list[[i]] <- vcv.phylo(phy=mytrees[[i]],mydata.mean.rs)
    } 
    
#Simulate trait evolution on each phylogeny
  shape.sim <- as.list(rep(NA,length(mytrees)))
    
    for (i in 1: length(mytrees)){
      shape.sim[[i]] <- sim.char(mytrees[[i]], varcov.list[[i]], nsim=1000, model="BM")
    }
    
#Combine simulations into one list
  simlist <- list.arrays.to.matrices(shape.sim)

#----------------------------------------------------------
#Compare observed and simulated data

#PCA on observed data
  mydata.mean.PCA <- prcomp(mydata.mean.rs)
#Select PC axes that account for 95% of the variation
  mydata.mean.PCaxes <- selectPCaxes.prcomp(mydata.mean.PCA, 0.956)
  #Calculate disparity as sum of variance
  mydata.mean.sumvar <- PCsumvar(mydata.mean.PCaxes) 
  
#PCA on simulated data
  simlist.PCA <- NULL
  
  for (i in 1:length(simlist)){
    simlist.PCA[[i]] <- prcomp(simlist[[i]])
  }

#Select PC axes that account for 95% of the variation
  simlist.PCaxes <- NULL
  
  for (i in 1:length(simlist.PCA)){
    simlist.PCaxes[[i]] <- selectPCaxes.prcomp(simlist.PCA[[i]], 0.956)
  }
  
#Calculate disparity as sum of variance
  simlist.sumvar <- NULL
  
  for (i in 1:length(simlist.PCA)){
    simlist.sumvar[[i]] <- PCsumvar(simlist.PCaxes[[i]])
  }

#Compare the observed and simulated values  
p.sumvar <- pvalue.dist(simlist.sumvar, mydata.mean.sumvar)

#Histogram comparison
  sumvar.hist <- hist(simlist.sumvar, xlab="Sum of Variance", main=NULL, las=1,
                      cex.lab=1.2)
    arrow.to.x.point(sumvar.hist, mydata.mean.sumvar, fraction.of.yaxis=50, line.fraction.of.yaxis=4,
                    height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)