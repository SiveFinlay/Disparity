#22/07/2014

#New simulation analyses after Dean Adams' email advice;
  #My previous results had issues with scale: the simulated values were not of a comparable scale to my original data

#He gave two options: one based on ratios and the other on "back transformation" - this script is for back transformation

# Find a way to ‘back scale’ the simulated datasets to the original Procrustes units and compare each group’s disparity separately.
    # This latter approach requires some thought to implement correctly. It seems to me that since you are assuming BM, that is analogous to
    # simulating error in landmark locations for each landmark. So what you could do is take the vector of simulated values for each species,
    # add the overall reference specimen (overall mean from GPA) to these. Now you have simulated BM deviations from the overall mean form.
    # Then do GPA of the simulated dataset and calculate disparity. I believe this would put things back into Procrustes-comparable units.
    

library(ape)
library(geiger)
library(geomorph)


source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PValueFunction_FromDistribution.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r" )

#SkDors
#1) Phylogenies
   setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
   mytrees <- read.tree("SkDors_tenrec+gmole_101trees.phy")

#2) Data
   setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skdors")

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

#PCA of mean shape values for each species
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

#Ratios of observed tenrec/gmole disparity
  obs.sumvar.ratio <- tenrec.sumvar/gmole.sumvar
  obs.prodvar.ratio <- tenrec.prodvar/gmole.prodvar
  obs.sumrange.ratio <- tenrec.sumrange/gmole.sumrange
  obs.prodrange.ratio <- tenrec.prodrange/gmole.prodrange

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
      shape.sim[[i]] <- sim.char(mytrees[[i]], varcov[[i]],nsim=5,model="BM")         #scale up the simulations later
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

#Calculate disparity for each of these sets of Procrustes-superimposed simulated data
  sim.PCA <- NULL
    for (i in 1:length(sim.GPA)){
      sim.PCA[[i]] <- plotTangentSpace(sim.GPA[[i]]$coords)
    }

#Select PC axes
  sim.PC95axes <- NULL
    for (i in 1:length(sim.PCA)){
      sim.PC95axes[[i]] <- selectPCaxes(sim.PCA[[i]], 0.956, tax.sorted$Binomial)
    }
    #BREAKS AT NUMBER 167: FIRST PC axis is above the threshold so I need to fix that function