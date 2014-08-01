#04/07/2014
  #Updated on 15/07/14 after reading through Luke Mahler's Continuous_tutorial script (under R_work/tutorials)

#Script to test my simulations code on linear measurements
  #Data file is my mandible measurements
  #(created by taking the intact mandible measurements from Skulls_after_remeasuring_06_2013
  #and combining it with the mandible measurements from Skulls_FMNH_Sept2013)
  


library(ape)
library(geiger)
library(phytools)
library(qpcR)

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
    mydata.mean.rs[,i] <- as.character(mydata.mean$Mean[which(mydata.mean$Measure == (levels(mydata.mean$Measure)[i]))])
   }                              #without as.character then the code just ranks the values

#Log each of the variables
  mydata.mean.rs.log <- matrix(NA, length(binom), length(measures))
  
  for (i in 1:length(measures)){
    rownames(mydata.mean.rs.log) <- binom
    colnames(mydata.mean.rs.log) <- levels(mydata.mean$Measure)
    mydata.mean.rs.log[,i] <- log(as.numeric(mydata.mean.rs[,i]))
  }
############################################
#NB; I missed a possibly very important step: sort the trait data to be in the same order as the tip labels
  #(It won't make a difference if the code can recognise the same binomial species names in the data and the trees
  # but Mahler's tutorial code advocates that it's a good idea to sort data)
  
#However, there's an issue in that the tip data are in different orders in every tree
#So presumably I need to create differently sorted trait data matrices to match each phylogeny

mydata.sorted <- NULL

for (i in 1:length(mytrees)){
  mydata.sorted[[i]] <- mydata.mean.rs.log[mytrees[[i]]$tip.label,]
}


#-------------------------------------------------------
#Variance covariance matrix for the trait data
  #NB: advice from Adam Algar at BES Macro Nottingham: use ic.sigma
      #but that's been deprecated -> use vcv.phylo instead (previously I used vcv)
      #Here using vcv and vcv.phylo give the same results

#One phylogeny first
 one.tree <- mytrees[[1]]
 
 varcov.phylo <- vcv.phylo(one.tree, mydata.mean.rs.log)
  varcov <- vcv(one.tree, mydata.mean.rs.log)
  varcov.phylo==varcov #All true
#-----------------------------------------------------------
#Separate variance covariance matrix for each of the phylogenies

#Separate variance covariance matrix of the logged average trait values for each of the phylogenies
  #use the list of sorted trait matrices
  varcov.list <- as.list(rep(NA,length(mytrees)))
    
    for(i in 1:length(mytrees)){
      varcov.list[[i]] <- vcv.phylo(phy=mytrees[[i]],mydata.sorted[[i]])
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

#PCA on observed data (log values)
  mydata.mean.PCA <- prcomp(mydata.mean.rs.log)
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
                    
#---------------------------------
#Observed disparity is significantly lower than the simulated values, 
  #even after sorting the data according to the tip labels in the phylogeny

#----------------------------------------------
#Trying out code from Mahler's Continuous_tutorial script

#Use the PCA results of the trait values
#Select PC1 (drop=FALSE maintains the same object class - stops it from turning a single column into a vector)

mands.pc <- (mydata.mean.PCA$x[,1, drop=FALSE])

#Working with one tree only

#Check whether the data and species in the tree are in the same order
(mytrees[[1]]$tip.label) == rownames(mands.pc)   #not true for all of them

#Re-order the data based on the order of the tip labels in one phylogeny only
    mands.pc.order <- mands.pc[one.tree$tip.label,, drop=FALSE]
    
    
#Re-scale variables to plot them to look for variation
	mands.lines <- (mands.pc-min(mands.pc)) / (max(mands.pc)-min(mands.pc))

#Plot the single tree in one half of the plotting window	
	plot(one.tree, x.lim=250,cex=.5,font=4,label.offset=.025,edge.width=3)  

#plot line segments corresponding to each trait.

	nspecies <- length(one.tree$tip.label)
	segments(rep(150,nspecies),1:nspecies,rep(150,nspecies)+(25*mands.lines), 1:nspecies,lwd=3)
	mtext("relative mandible length", at = 150, side = 1, line = 0, cex=.5,font=2)

#Try different models of continuous trait evolution	
  brown_mands <- fitContinuous(one.tree, mands.pc.order)	
	eb_mands <- fitContinuous(one.tree, mands.pc.order, model="EB") #using the default bounds instead of the tutorial example, 
                                                                                          #I don't understand how they're chosen
	ou_mands <- fitContinuous(one.tree, mands.pc.order, ,model="OU")

#Table to compare the models	
	model<-matrix(,3,4,dimnames = list(c("Brownian Motion", "Early Burst", "Ornstein-Uhlenbeck"),c("log likelihood", "AICc", "Delta AICc", "AICc Weights")))

# Now, let's put the likelihood scores and AICc values into the first two columns, calling them from the dataframes we made earlier:

	model[,1]<-c(brown_mands$opt$lnL, eb_mands$opt$lnL, ou_mands$opt$lnL)
	model[,2]<-c(brown_mands$opt$aicc, eb_mands$opt$aicc, ou_mands$opt$aicc)

# Add the delta AIC scores and the AIC weights 

	aic.all<-as.matrix(model[,2])
	scor.wts<-akaike.weights(aic.all)


	model[,3]<-scor.wts$deltaAIC
	model[,4]<-scor.wts$weights

  #NB: this method works but it's not particularly informative since there isn't an obvious difference in AIC weights
    #and also I've only tried to fo it on one tree; it would take longer to fit all models to all of the trees

