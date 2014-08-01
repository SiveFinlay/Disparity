#22/07/2014

#New simulation analyses after Dean Adams' email advice;
  #My previous results had issues with scale: the simulated values were not of a comparable scale to my original data

#His advice was to compare the observed differences between tenrec and gmole disparity to simulations under a BM model
  #So I'm not comparing tenrecs on their own to chance but rather the ratio of disp_tenrecs/disp_gmoles to the simulated values of that ratio
  
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

#------------------------------------------
#Disparity calculations for simulated data

#PCA
  shape.simPC <- calc.each.list(mylist=simlist.sorted, calculation=prcomp)

#Select the PC axes that account for approximately 95% of the variation
  shape.simPC95 <- NULL
    for (i in 1:length(shape.simPC)){
      shape.simPC95[[i]] <- selectPCaxes.prcomp(shape.simPC[[i]], 0.956)
    }

#Calculate group disparity for each simulation
#Rows that correspond to tenrecs
  sim.tenrec.rows <- which(tax.sorted$Family=="Tenrecidae")

  #Select just the tenrec rows from each simulated data set
  sim.tenrec <- NULL

    for (i in 1:length(shape.simPC95)){
      sim.tenrec[[i]] <- shape.simPC95[[i]][sim.tenrec.rows,]
    }

   sim.tenrec.sumvar <- unlist(calc.each.list(mylist=sim.tenrec, calculation=PCsumvar))
   sim.tenrec.prodvar <- unlist(calc.each.list(mylist=sim.tenrec, calculation=PCprodvar))
   
   sim.tenrec.sumrange <- unlist(calc.each.list(mylist=sim.tenrec, calculation=PCsumrange))
   sim.tenrec.prodrange <- unlist(calc.each.list(mylist=sim.tenrec, calculation=PCprodrange))
   
#Rows that correspond to golden moles
  sim.gmole.rows <- which(tax.sorted$Family=="Chrysochloridae")

  #Select just the gmole rows from each simulated data set
  sim.gmole <- NULL

    for (i in 1:length(shape.simPC95)){
      sim.gmole[[i]] <- shape.simPC95[[i]][sim.gmole.rows,]
    }

   sim.gmole.sumvar <- unlist(calc.each.list(mylist=sim.gmole, calculation=PCsumvar))
   sim.gmole.prodvar <- unlist(calc.each.list(mylist=sim.gmole, calculation=PCprodvar))

   sim.gmole.sumrange <- unlist(calc.each.list(mylist=sim.gmole, calculation=PCsumrange))
   sim.gmole.prodrange <- unlist(calc.each.list(mylist=sim.gmole, calculation=PCprodrange))

#Ratios of simulated tenrec/gmole disparity
  sim.sumvar.ratio <- sim.tenrec.sumvar/sim.gmole.sumvar
  sim.prodvar.ratio <- sim.tenrec.prodvar/sim.gmole.prodvar
  sim.sumrange.ratio <- sim.tenrec.sumrange/sim.gmole.sumrange
  sim.prodrange.ratio <- sim.tenrec.prodrange/sim.gmole.prodrange
  
#--------------------------------------------------------
#Compare the observed ratios to the simulated ratios


  sumvar.p <- pvalue.dist(distribution=sim.sumvar.ratio, obs.val=obs.sumvar.ratio)
  prodvar.p <- pvalue.dist(distribution=sim.prodvar.ratio, obs.val=obs.prodvar.ratio)
  sumrange.p <- pvalue.dist(distribution=sim.sumrange.ratio, obs.val=obs.sumrange.ratio)
  prodrange.p <- pvalue.dist(distribution=sim.prodrange.ratio, obs.val=obs.prodrange.ratio)
