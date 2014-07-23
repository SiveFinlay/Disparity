#22/07/2014

#New simulation analyses after Dean Adams' email advice;
  #My previous results had issues with scale: the simulated values were not of a comparable scale to my original data

#He gave two options: one based on ratios and the other on "back transformation" - this script is for back transformation

# Find a way to ‘back scale’ the simulated datasets to the original Procrustes units and compare each group’s disparity separately.
    # This latter approach requires some thought to implement correctly. It seems to me that since you are assuming BM, that is analogous to
    # simulating error in landmark locations for each landmark. So what you could do is take the vector of simulated values for each species,
    # add the overall reference specimen (overall mean from GPA) to these. Now you have simulated BM deviations from the overall mean form.
    # Then do GPA of the simulated dataset and calculate disparity. I believe this would put things back into Procrustes-comparable units.

#This method seems to work but I need to up-scale it to more simulations
  #However first I need to sort out the problem of dealing with different versions of geiger
      
#-------------------------------------------------------
library(ape)
library(geiger)   
#NB: The new version of geiger (geiger_2.0.3.)creates different outputs for the vcv.phylo function
      #it produces vcv matrices with species as the rows and columns, not traits
      #I have an older geiger version (geiger_1.99-3.1) installed on my laptop but the new version is the alien
library(geomorph)


source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PValueFunction_FromDistribution.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r" )

#On cyberman (alien as remote server)
  #source("~/Disparity/functions/Disparity_general_functions.r")
  #source("~/Disparity/functions/DisparityFunctions_Variance_Range.r")
  #source("~/Disparity/functions/PvalueFunction_FromDistribution.r")





#SkDors
#1) Phylogenies
   setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
   #setwd("~/Disparity/output/phylogenies")
   mytrees <- read.tree("SkDors_tenrec+gmole_101trees.phy")

#2) Data
   setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skdors")
   #setwd("~/Disparity/output/shape_data/skdors")
   
   
  #2a) All tenrecs and golden moles
      #shape coordinates
      sps.mean <- dget(file="SkDors_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      tax <- read.table("SkDors_tenrec+gmole_sps.mean_taxonomy.txt")


#Don't need to prune the phylogenies or change the data because it's just golden moles and tenrecs

#################################
#VCV matrix of superimposed Procrustes coordinates

#Convert the shape coordinates into a 2D array
  twoDshape <- two.d.array(sps.mean$meanshape)

#Add the species as rownames
  rownames(twoDshape) <- sps.mean$Binom

#-------------------------------------------------
#Calculate observed disparity

#PCA of mean shape values for each species     (NB: L.Harmon said to do a phylogenetic PCA)
  obsPC <- prcomp(twoDshape)

#select the PC axes which correspond to 95% of the variation
  obsPC95 <- selectPCaxes.prcomp(obsPC, 0.956)


#Observed disparity for tenrecs
  tenrecPC <- obsPC95[which(tax$Family=="Tenrecidae"),]

  tenrec.sumvar <- PCsumvar(tenrecPC)
  tenrec.prodvar <- PCprodvar(tenrecPC)
  tenrec.sumrange <- PCsumrange(tenrecPC)
  tenrec.prodrange <- PCprodrange(tenrecPC)

#Observed disparity for golden moles
  gmolePC <- obsPC95[which(tax$Family=="Chrysochloridae"),]

  gmole.sumvar <- PCsumvar(gmolePC)
  gmole.prodvar <- PCprodvar(gmolePC)
  gmole.sumrange <- PCsumrange(gmolePC)
  gmole.prodrange <- PCprodrange(gmolePC)

#---------------------------------------------------------------
#Simulate data from VCV matrices

#NB: Sort the data frames so that the species are in the same order as the tip labels on the phylogenies
  #(This may not be important but it's worth doing just in case, cf. Mahler's Continuous_tutorial script)
  twoDshape.sorted <- NULL

  for (i in 1:length(mytrees)){
    twoDshape.sorted[[i]] <- twoDshape[mytrees[[i]]$tip.label,]
  }

#Use the sorted data that corresponds with each tree

#Separate variance covariance matrix of the shape data for each of the phylogenies
  varcov <- as.list(rep(NA,length(mytrees)))
    for(i in 1:length(mytrees)){
      varcov[[i]] <- vcv.phylo(phy=mytrees[[i]],twoDshape.sorted[[i]])
    }

#Simulate shape evolution on each phylogeny
  shape.sim <- as.list(rep(NA,length(mytrees)))
    for (i in 1: length(mytrees)){
      shape.sim[[i]] <- sim.char(mytrees[[i]], varcov[[i]],nsim=50,model="BM")         #scale up the simulations later
    }

#Combine simulations into one list
simlist <- list.arrays.to.matrices(shape.sim)

#Sort all of the matrices into one order

#Sort the taxonomic data alphabetically by species name
  tax.sorted <- tax[order(tax$Binomial),]

#Sort all of the simulated matrices into the same alphabetical order

    simlist.sorted <- NULL

    for (i in 1:length(simlist)){
      simlist.sorted[[i]] <- simlist[[i]][order(rownames(simlist[[i]])),]
    }

#-------------------------------------------
#Convert the simulated data back into arrays of landmark data

  sim.array <- NULL
    for (i in 1:length(simlist.sorted)){
      sim.array[[i]] <- arrayspecs(simlist.sorted[[i]],dim(sps.mean$meanshape)[1],2)
    }

#Add the overall reference specimen to each species

#Find the observed mean shape
  mean.shape <- mshape(sps.mean$meanshape)

#Add the mean shape to the simulated values for each species

  sim.add.mshape <- rep(list(array(NA, dim(sim.array[[i]]))), length(sim.array))
  
    for (i in 1:length(sim.array)){
      for (j in 1:dim(sim.array[[i]])[3]){
        sim.add.mshape[[i]][,,j] <- sim.array[[i]][,,j] + mean.shape
      }
    }
    
#-------------------------------------
#GPA of each simulated data set
  sim.GPA <- NULL

    for (i in 1:length(sim.add.mshape)){
      sim.GPA[[i]] <- gpagen(sim.add.mshape[[i]], ProcD=TRUE)
    }

#PCA for each of these sets of Procrustes-superimposed simulated data  
  #(NB: This takes ages to run because it creates a PCA graph for each one; maybe it would be better to use prcomp?)
  sim.PCA <- NULL
    for (i in 1:length(sim.GPA)){
      sim.PCA[[i]] <- plotTangentSpace(sim.GPA[[i]]$coords)
    }

#Select PC axes
  sim.PC95axes <- NULL
    for (i in 1:length(sim.PCA)){
      sim.PC95axes[[i]] <- selectPCaxes(sim.PCA[[i]], 0.956, tax.sorted$Binomial)
    }

#Calculate disparity for each simulated data set
  #Tenrecs
  tenrec.rows <- which(tax.sorted$Family == "Tenrecidae") 
  
  #tenrec rows of the PC matrices
    tenrec.sim.PC <- NULL
      for (i in 1:length(sim.PC95axes)){
        tenrec.sim.PC[[i]] <- sim.PC95axes[[i]][tenrec.rows,]
      }
   
  #Calculate tenrec disparity   
         
  sim.tenrec.sumvar <- unlist(calc.each.list(mylist=tenrec.sim.PC, calculation=PCsumvar))
  sim.tenrec.prodvar <- unlist(calc.each.list(mylist=tenrec.sim.PC, calculation=PCprodvar))
   
  sim.tenrec.sumrange <- unlist(calc.each.list(mylist=tenrec.sim.PC, calculation=PCsumrange))
  sim.tenrec.prodrange <- unlist(calc.each.list(mylist=tenrec.sim.PC, calculation=PCprodrange))


  #Golden moles
    gmole.rows <- which(tax.sorted$Family == "Chrysochloridae") 
  
  #golden mole rows of the PC matrices
    gmole.sim.PC <- NULL
      for (i in 1:length(sim.PC95axes)){
        gmole.sim.PC[[i]] <- sim.PC95axes[[i]][gmole.rows,]
      }
   
  #Calculate gmoledisparity   
         
  sim.gmole.sumvar <- unlist(calc.each.list(mylist=gmole.sim.PC, calculation=PCsumvar))
  sim.gmole.prodvar <- unlist(calc.each.list(mylist=gmole.sim.PC, calculation=PCprodvar))
   
  sim.gmole.sumrange <- unlist(calc.each.list(mylist=gmole.sim.PC, calculation=PCsumrange))
  sim.gmole.prodrange <- unlist(calc.each.list(mylist=gmole.sim.PC, calculation=PCprodrange))

#----------------------------------------------
#Compare observed and simulated disparity
  #Tenrecs   
    tenrec.sumvar.p <- pvalue.dist(distribution=sim.tenrec.sumvar, obs.val=tenrec.sumvar)
    tenrec.prodvar.p <- pvalue.dist(distribution=sim.tenrec.prodvar, obs.val=tenrec.prodvar) 
    tenrec.sumrange.p <- pvalue.dist(distribution=sim.tenrec.sumrange, obs.val=tenrec.sumrange) 
    tenrec.prodrange.p <- pvalue.dist(distribution=sim.tenrec.prodrange, obs.val=tenrec.prodrange) 
    
  #Golden moles   
    gmole.sumvar.p <- pvalue.dist(distribution=sim.gmole.sumvar, obs.val=gmole.sumvar)
    gmole.prodvar.p <- pvalue.dist(distribution=sim.gmole.prodvar, obs.val=gmole.prodvar) 
    gmole.sumrange.p <- pvalue.dist(distribution=sim.gmole.sumrange, obs.val=gmole.sumrange) 
    gmole.prodrange.p <- pvalue.dist(distribution=sim.gmole.prodrange, obs.val=gmole.prodrange)   