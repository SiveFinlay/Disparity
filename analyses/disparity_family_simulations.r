#16/07/2014 (date of most recent modification)
#General script for simulating shape evolution across phylogenies and calculating disparity
    #Modified the script to make it work for running the analyses through the terminal connection to cyberman
    #Need to re-run the analysis to check whether sorting the data makes any difference
    
#Steps
  #1) Read in phylogenies, shape data and taxonomy
  #2) Choose which family (tenrecs or golden moles) 
  #3) Prune phylogenies
  #Extra step: sort the data to match the order of the tip labels in each phylogeny
  #4) Shape simulation across phylogenies
  #5) PCA analysis of each simulation
  #6) Calculate disparity for each simulation and observed data
  #7) Compare observed and simulated disparity
  #8) Create output files
        #Table of disparity comparisons (observed vs. simulated)
        #Histograms of disparity comparisons; sum of variance
                                            # product of variance
                                            # sum of ranges 
                                            # product of ranges 
                                            # ZelditchMD 
  
#########################################################
library(ape)
library(geiger)
library(geomorph)

#-------------------------------------------------------------------------------------
#First option: working directories

#Run on my computer
	#source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
	#source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PValueFunction_FromDistribution.r")
  #source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r" )

#On the alien: save everything onto the USB
  #source("E:/Disparity/functions/Disparity_general_functions.r")
  #source("E:/Disparity/functions/DisparityFunctions_Variance_Range.r")
  #source("E:/Disparity/functions/PValueFunction_FromDistribution.r")
  
#On cyberman (alien as remote server)
  source("~/Disparity/functions/Disparity_general_functions.r")
  source("~/Disparity/functions/DisparityFunctions_Variance_Range.r")
  source("~/Disparity/functions/PvalueFunction_FromDistribution.r")

######################################################
#1) READ IN DATA
######################################################

#SkDors
#1) Phylogenies
   #setwd("~/Disparity/output/phylogenies")
   #mytrees <- read.tree("SkDors_tenrec+gmole_101trees.phy")

#2) Data
   #setwd("~/Disparity/output/shape_data/skdors")
  
  #2a) All tenrecs and golden moles
      #shape coordinates
      #sps.mean <- dget(file="SkDors_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      #tax <- read.table("SkDors_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
      #shape coordinates
      #sps.mean <- dget(file="SkDors_nonmic_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      #tax <- read.table("SkDors_nonmic_tenrec+gmole_sps.mean_taxonomy.txt")
#------------------------------------------------------
#SkLat

#1) Phylogenies
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
setwd("~/Disparity/output/phylogenies")
     mytrees <- read.tree("SkLat_tenrec+gmole_101trees.phy")
     
#2) Data
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/sklat")
setwd("~/Disparity/output/shape_data/sklat")
  
  #2a) All tenrecs and golden moles
      #shape coordinates
      sps.mean <- dget(file="SkLat_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      tax <- read.table("SkLat_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
    #sps.mean <- dget(file="SkLat_nonmic_tenrec+gmole_sps.mean.txt")
    #taxonomic information
    #tax <- read.table("SkLat_nonmic_tenrec+gmole_sps.mean_taxonomy.txt")
#------------------------------------------------------
#SkVent
#1) Phylogenies
#setwd("~/Disparity/output/phylogenies")
    #mytrees <- read.tree("SkVent_tenrec+gmole_101trees.phy")

#2) Data
#setwd("~/Disparity/output/shape_data/skvent")
  
  #2a) All tenrecs and golden moles
      #shape coordinates
      #sps.mean <- dget(file="SkVent_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      #tax <- read.table("SkVent_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
    #sps.mean <- dget("SkVent_nonmic_tenrec+gmole_sps.mean.txt")
    #taxonomic information
    #tax <- read.table("SkVent_nonmic_tenrec+gmole_sps.mean_taxonomy.txt")
#------------------------------------------------------
#Mandibles
#1) Phylogenies
#setwd("~/Disparity/output/phylogenies")
    #mytrees <- read.tree("Mands_tenrec+gmole_101trees.phy")

#2) Data
#setwd("~/Disparity/output/shape_data/mands")
  #2a) All tenrecs and golden moles
      #shape coordinates
      #sps.mean <- dget(file="Mands_tenrec+gmole_sps.mean.txt")
      #taxonomic information
      #tax <- read.table("Mands_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
    #sps.mean <- dget(file="Mands_nonmic_tenrec+gmole_sps.mean.txt")
    #taxonomic information
    #tax <- read.table("Mands_nonmic_tenrec+gmole_sps.mean_taxonomy.txt")

#################################################
#2) CHOOSE WHICH FAMILY 
#################################################
#Golden moles
  #fam <- "Chrysochloridae"
#Tenrecs
  fam <- "Tenrecidae"

  sps.tax <- tax$Binomial[which(tax$Family == fam)]

#find the ID numbers for the species of interest
  ID.sps <- matching.id(sps.tax, sps.mean$Binom)

#select those species from the overall data
  mysps.mean <- select.from.list(sps.mean, ID.sps)
#drop unused levels
  mysps.mean <- droplevels.from.list(mysps.mean)

##################################################
#3) PRUNE THE PHYLOGENIES
##################################################
#Prune the trees to include that family's taxa only

  TreeOnly <- tree.only(mytrees, sps.tax)

#prune the trees so that they only include the species which are in the species data
  sps.trees <- remove.missing.species.tree(mytrees, TreeOnly)

###################################################
#4) SHAPE SIMULATION
###################################################

#Convert the shape coordinates into a 2D array
  twoDshape <- two.d.array(mysps.mean$meanshape)

#Add the species as rownames
  rownames(twoDshape) <- mysps.mean$Binom

#NB: Sort the data frames so that the species are in the same order as the tip labels on the phylogenies
  #(cf. Mahler's Continuous_tutorial script
  twoDshape.sorted <- NULL

  for (i in 1:length(sps.trees)){
    twoDshape.sorted[[i]] <- twoDshape[sps.trees[[i]]$tip.label,]
  }

#Use the sorted data that corresponds with each tree

#Separate variance covariance matrix of the shape data for each of the phylogenies
  varcov <- as.list(rep(NA,length(sps.trees)))
    for(i in 1:length(sps.trees)){
      varcov[[i]] <- vcv.phylo(phy=sps.trees[[i]],twoDshape.sorted[[i]])
    }   

#simulate shape evolution on each phylogeny
  shape.sim <- as.list(rep(NA,length(sps.trees)))
    for (i in 1: length(sps.trees)){
      shape.sim[[i]] <- sim.char(sps.trees[[i]], varcov[[i]],nsim=1000,model="BM")
    }
    
#Combine simulations into one list
simlist <- list.arrays.to.matrices(shape.sim)  
######################################################
#5) PCA ANALYSIS
######################################################

#a) Simulated data
  shape.simPC <- calc.each.list(mylist=simlist, calculation=prcomp)

#Select the PC axes that account for approximately 95% of the variation
  shape.simPC95 <- NULL
    for (i in 1:length(shape.simPC)){
      shape.simPC95[[i]] <- selectPCaxes.prcomp(shape.simPC[[i]], 0.956)
    }
  
#b) Observed data
  #Do a principal components analysis on the family's (tenrec or golden mole)shape values only
    #i.e. don't use the PC axes from the global principal components analysis of all the species together
    #Makes sense because simulations are based on the shape coordinates for one family only

#PCA of mean shape values for each species
  obsPC <- prcomp(twoDshape)

#select the PC axes which correspond to 95% of the variation
  obsPC95 <- selectPCaxes.prcomp(obsPC, 0.956)

##########################################################
#6) CALCULATE DISPARITY FOR SIMULATIONS AND OBSERVED DATA
##########################################################
#Simulated data
#a) Disparity based on PC axes
  #Simulated data
    #Variance measures
    sumvar <- calc.each.list(mylist=shape.simPC95, calculation=PCsumvar)
    prodvar <- calc.each.list(mylist=shape.simPC95, calculation=PCprodvar)
    
      #Matrix of the sum and product of variance for each simulation
      simPC.var <- as.data.frame(matrix(NA,nrow=length(sumvar),ncol=2))
        colnames(simPC.var) <- c("SumVar","ProdVar")
        simPC.var[,1] <- unlist(sumvar)
        simPC.var[,2] <- unlist(prodvar)
  

    #Range measures
    sumrange <- calc.each.list(mylist=shape.simPC95, calculation=PCsumrange)
    prodrange <- calc.each.list(mylist=shape.simPC95, calculation=PCprodrange)
      
      #Matrix of the sum and product of range for each simulation
      simPC.range <- as.data.frame(matrix(NA,nrow=length(sumrange),ncol=2))
        colnames(simPC.range) <- c("SumRange","ProdRange")
        simPC.range[,1] <- unlist(sumrange)
        simPC.range[,2] <- unlist(prodrange)

  #Observed data
    #Variance
    obs.sumvar <- PCsumvar(obsPC95)
    obs.prodvar <- PCprodvar(obsPC95)
  
    #Range
    obs.sumrange <- PCsumrange(obsPC95)
    obs.prodrange <- PCprodrange(obsPC95)

#b) Disparity based on interlandmark distances
  #Simulated data
    #convert the simulated shape matrices into three dimensional arrays
      simlist.arrays <- NULL
        for (i in 1:length(simlist)){
          simlist.arrays[[i]] <- arrayspecs(A=simlist[[i]], p=((dim(simlist[[i]])[2])/2), k=2)
        }                                            #(Old version of geomorph had a byLAND=FALSE option)

    #calculate the ild.distances for each simulation: compare each species to the overall mean shape of all species
      simlist.ild <- NULL
        for (i in 1:length(simlist.arrays)){
          simlist.ild[[i]]<-dist.to.ref(simlist.arrays[[i]])
        }

    #calculate disparity as the sum of squared distances (Zeldich 2012)
      sim.md <- NULL
        for (i in 1:length(simlist.ild)){
          sim.md[[i]] <- ZelditchMD(simlist.ild[[i]])
        }

  #Observed data
    #calculate the ild.distance for each species
      obs.ild <- dist.to.ref(mysps.mean$meanshape) 
    
    #calculate disparity
      obs.md <- ZelditchMD(obs.ild)


######################################################### 
#7)COMPARE OBSERVED AND SIMULATED DISPARITY
#########################################################

#Compare observed disparity to the distribution of simulated values
  # (histograms in the output section below)
  sumvar.p <- pvalue.dist(distribution=simPC.var$SumVar, obs.val=obs.sumvar)
  prodvar.p <- pvalue.dist(distribution=simPC.var$ProdVar, obs.val=obs.prodvar)
  sumrange.p <- pvalue.dist(distribution=simPC.range$SumRange, obs.val=obs.sumrange)
  prodrange.p <- pvalue.dist(distribution=simPC.range$ProdRange, obs.val=obs.prodrange)
  
  md.p <- pvalue.dist(distribution=sim.md, obs.val=obs.md)


#Create a table to compare the disparity measures
  disp <- as.data.frame(matrix(NA, nrow=5, ncol=5))
  rownames(disp) <- c("SumVar","ProdVar","SumRange","ProdRange","ZelditchMD")
  colnames(disp) <- c("Observed","Sim.min","Sim.max", "Sdev.sim","p.value")

    disp[1,1] <- obs.sumvar
    disp[1,2] <- range(simPC.var$SumVar)[1]
    disp[1,3] <- range(simPC.var$SumVar)[2]
    disp[1,4] <- sd(simPC.var$SumVar)
    disp[1,5] <- sumvar.p

    disp[2,1] <- obs.prodvar
    disp[2,2] <- range(simPC.var$ProdVar)[1]
    disp[2,3] <- range(simPC.var$ProdVar)[2]
    disp[2,4] <- sd(simPC.var$ProdVar)
    disp[2,5] <- prodvar.p

    disp[3,1] <- obs.sumrange
    disp[3,2] <- range(simPC.range$SumRange)[1]
    disp[3,3] <- range(simPC.range$SumRange)[2]
    disp[3,4] <- sd(simPC.range$SumRange)
    disp[3,5] <- sumrange.p

    disp[4,1] <- obs.prodrange
    disp[4,2] <- range(simPC.range$ProdRange)[1]
    disp[4,3] <- range(simPC.range$ProdRange)[2]
    disp[4,4] <- sd(simPC.range$ProdRange)
    disp[4,5] <- prodrange.p
    
    disp[5,1] <- obs.md
    disp[5,2] <- range(sim.md)[1]
    disp[5,3] <- range(sim.md)[2]
    disp[5,4] <- sd(sim.md)
    disp[5,5] <- md.p

#######################################
#8) CREATE THE OUTPUT FILES
#######################################
#Set the correct working directory for each data set 

#SkDors
  #setwd("~/Disparity/output/shape_simulations/skdors")

#SkLat
  setwd("~/Disparity/output/shape_simulations/sklat")
  
#SkVent
  #setwd("~/Disparity/output/shape_simulations/skvent")
  
#Mands
  #setwd("~/Disparity/output/shape_simulations/mands")

#*******************************************************************************
# Table of disparity comparisons
#*******************************************************************************
#SkDors
  #Tenrecs
  #write.table(file="skdors_trc+gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #Golden moles
  #write.table(file="skdors_trc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

  #Non-microgale tenrecs
  #write.table(file="skdors_nonmictrc+gmole_nonmic_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#-------------------------------------------------------------------------------  
#SkLat
  #Tenrecs
  #write.table(file="sklat_trc+gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #Golden moles
  #write.table(file="sklat_trc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

  #Non-Microgale tenrecs
  write.table(file="sklat_nonmictrc+gmole_nonmic_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#-------------------------------------------------------------------------------
#SkVent
  #Tenrecs
  #write.table(file="skvent_trc+gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #Golden moles
  #write.table(file="skvent_trc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

  #Non-Microgale tenrecs
  #write.table(file="skvent_nonmictrc+gmole_nonmic_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#-------------------------------------------------------------------------------
#Mands
  #Tenrecs
  #write.table(file="mands_trc+gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #Golden moles
  #write.table(file="mands_trc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

  #Non-Microgale tenrecs
  #write.table(file="mands_nonmictrc+gmole_nonmic_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#********************************************************
#Histogram plots of the simulated disparity values
  #arrows point to the observed disparity values for comparison
  #breaks in the x axes because the limits are extended to include the observed values
  #saved separate plots for each comparison because otherwise they're too small to see
#*******************************************************
#1/5) Sum of variance 

#SkDors
  #Tenrecs
  #pdf(file="skdors_trc+gmole_tenrec_sumvariance.pdf")
  #Golden moles
  #pdf(file="skdors_trc+gmole_gmole_sumvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skdors_nonmictrc+gmole_nonmic_tenrec_sumvariance.pdf")
#-------------------------------------------------------------------
#SkLat
  #Tenrecs
  #pdf(file="sklat_trc+gmole_tenrec_sumvariance.pdf")
  #Golden moles
  #pdf(file="sklat_trc+gmole_gmole_sumvariance.pdf")
  
  #Non-microgale tenrecs
  pdf(file="sklat_nonmictrc+gmole_nonmic_tenrec_sumvariance.pdf")
#-------------------------------------------------------------------
#SkVent
  #Tenrecs
  #pdf(file="skvent_trc+gmole_tenrec_sumvariance.pdf")
  #Golden moles
  #pdf(file="skvent_trc+gmole_gmole_sumvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skvent_nonmictrc+gmole_nonmic_tenrec_sumvariance.pdf")

#--------------------------------------------------------------------
#Mands
  #Tenrecs
  #pdf(file="mands_trc+gmole_tenrec_sumvariance.pdf")
  #Golden moles
  #pdf(file="mands_trc+gmole_gmole_sumvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="mands_nonmictrc+gmole_nonmic_tenrec_sumvariance.pdf")
#------------------------------------------------------------
par(mgp=c(3.2,0.5,0))
par(mar=c(5,7,4,2)+0.1)                   
 sumvar.hist <- hist(simPC.var$SumVar, xlab="Sum of Variance", main=NULL, las=1, 
                  xlim=c(obs.sumvar, max(simPC.var$SumVar)), cex.lab=1.2)
    arrow.to.x.point(sumvar.hist, obs.sumvar, fraction.of.yaxis=50, line.fraction.of.yaxis=4, 
                     height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)
                     
  dev.off()  
  
#**************************************************************
#2/5) Product of variance 

#SkDors
  #Tenrecs
  #pdf(file="skdors_trc+gmole_tenrec_prodvariance.pdf")
  #Golden moles
  #pdf(file="skdors_trc+gmole_gmole_prodvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skdors_nonmictrc+gmole_nonmic_tenrec_prodvariance.pdf")
#-------------------------------------------------------------------
#SkLat
  #Tenrecs
  #pdf(file="sklat_trc+gmole_tenrec_prodvariance.pdf")
  #Golden moles
  #pdf(file="sklat_trc+gmole_gmole_prodvariance.pdf")
  
  #Non-microgale tenrecs
  pdf(file="sklat_nonmictrc+gmole_nonmic_tenrec_prodvariance.pdf")
#-------------------------------------------------------------------
#SkVent
  #Tenrecs
  #pdf(file="skvent_trc+gmole_tenrec_prodvariance.pdf")
  #Golden moles
  #pdf(file="skvent_trc+gmole_gmole_prodvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skvent_nonmictrc+gmole_nonmic_tenrec_prodvariance.pdf")

#--------------------------------------------------------------------
#Mands
  #Tenrecs
  #pdf(file="mands_trc+gmole_tenrec_prodvariance.pdf")
  #Golden moles
  #pdf(file="mands_trc+gmole_gmole_prodvariance.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="mands_nonmictrc+gmole_nonmic_tenrec_prodvariance.pdf")
#--------------------------------------------------------------------  
par(mgp=c(3.2,0.5,0))
par(mar=c(5,7,4,2)+0.1)                   
  prodvar.hist <- hist(simPC.var$ProdVar, xlab="Product of Variance", main=NULL, las=1, 
                  xlim=c(obs.prodvar, max(simPC.var$ProdVar)), cex.lab=1.2)
    arrow.to.x.point(prodvar.hist, obs.prodvar, fraction.of.yaxis=50, line.fraction.of.yaxis=4, 
                     height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)
                     
  dev.off()
#***************************************************************
#3/5) Sum of ranges 

#SkDors
  #Tenrecs
  #pdf(file="skdors_trc+gmole_tenrec_sumrange.pdf")
  #Golden moles
  #pdf(file="skdors_trc+gmole_gmole_sumrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skdors_nonmictrc+gmole_nonmic_tenrec_sumrange.pdf")
#-------------------------------------------------------------------
#SkLat
  #Tenrecs
  #pdf(file="sklat_trc+gmole_tenrec_sumrange.pdf")
  #Golden moles
  #pdf(file="sklat_trc+gmole_gmole_sumrange.pdf")
  
  #Non-microgale tenrecs
  pdf(file="sklat_nonmictrc+gmole_nonmic_tenrec_sumrange.pdf")
#-------------------------------------------------------------------
#SkVent
  #Tenrecs
  #pdf(file="skvent_trc+gmole_tenrec_sumrange.pdf")
  #Golden moles
  #pdf(file="skvent_trc+gmole_gmole_sumrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skvent_nonmictrc+gmole_nonmic_tenrec_sumrange.pdf")

#--------------------------------------------------------------------
#Mands
  #Tenrecs
  #pdf(file="mands_trc+gmole_tenrec_sumrange.pdf")
  #Golden moles
  #pdf(file="mands_trc+gmole_gmole_sumrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="mands_nonmictrc+gmole_nonmic_tenrec_sumrange.pdf")
#--------------------------------------------------------------------
   
par(mgp=c(3.2,0.5,0))
par(mar=c(5,7,4,2)+0.1)   
  sumrange.hist <- hist(simPC.range$SumRange, xlab="Sum of Ranges", main=NULL, las=1, 
                  xlim=c(obs.sumrange, max(simPC.range$SumRange)), cex.lab=1.2)
    arrow.to.x.point(sumrange.hist, obs.sumrange, fraction.of.yaxis=50, line.fraction.of.yaxis=4,
                    height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)

dev.off()

#***************************************************************
#4/5) Product of ranges 

#SkDors
  #Tenrecs
  #pdf(file="skdors_trc+gmole_tenrec_prodrange.pdf")
  #Golden moles
  #pdf(file="skdors_trc+gmole_gmole_prodrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skdors_nonmictrc+gmole_nonmic_tenrec_prodrange.pdf")
#-------------------------------------------------------------------
#SkLat
  #Tenrecs
  #pdf(file="sklat_trc+gmole_tenrec_prodrange.pdf")
  #Golden moles
  #pdf(file="sklat_trc+gmole_gmole_prodrange.pdf")
  
  #Non-microgale tenrecs
  pdf(file="sklat_nonmictrc+gmole_nonmic_tenrec_prodrange.pdf")
#-------------------------------------------------------------------
#SkVent
  #Tenrecs
  #pdf(file="skvent_trc+gmole_tenrec_prodrange.pdf")
  #Golden moles
  #pdf(file="skvent_trc+gmole_gmole_prodrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skvent_nonmictrc+gmole_nonmic_tenrec_prodrange.pdf")

#--------------------------------------------------------------------
#Mands
  #Tenrecs
  #pdf(file="mands_trc+gmole_tenrec_prodrange.pdf")
  #Golden moles
  #pdf(file="mands_trc+gmole_gmole_prodrange.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="mands_nonmictrc+gmole_nonmic_tenrec_prodrange.pdf")
#-------------------------------------------------------------------------
par(mgp=c(3.2,0.5,0))
par(mar=c(5,7,4,2)+0.1)                      
  prodrange.hist <- hist(simPC.range$ProdRange, xlab="Product of Ranges", main=NULL, las=1,
                    xlim=c(obs.prodrange, max(simPC.range$ProdRange)), cex.lab=1.2)
    arrow.to.x.point(prodrange.hist, obs.prodrange, fraction.of.yaxis=50, line.fraction.of.yaxis=4,
                    height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)
    
  dev.off()
#***************************************************************
#5/5) ZelditchMD (sum of squared distances)

#SkDors
  #Tenrecs
  #pdf(file="skdors_trc+gmole_tenrec_ZelditchMD.pdf")
  #Golden moles
  #pdf(file="skdors_trc+gmole_gmole_ZelditchMD.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skdors_nonmictrc+gmole_nonmic_tenrec_ZelditchMD.pdf")
  
#-------------------------------------------------------------------
#SkLat
  #Tenrecs
  #pdf(file="sklat_trc+gmole_tenrec_ZelditchMD.pdf")
  #Golden moles
  #pdf(file="sklat_trc+gmole_gmole_ZelditchMD.pdf")
  
  #Non-microgale tenrecs
  pdf(file="sklat_nonmictrc+gmole_nonmic_tenrec_ZelditchMD.pdf")
#-------------------------------------------------------------------
#SkVent
  #Tenrecs
  #pdf(file="skvent_trc+gmole_tenrec_ZelditchMD.pdf")
  #Golden moles
  #pdf(file="skvent_trc+gmole_gmole_ZelditchMD.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="skvent_nonmictrc+gmole_nonmic_tenrec_ZelditchMD.pdf")
#--------------------------------------------------------------------
#Mands
  #Tenrecs
  #pdf(file="mands_trc+gmole_tenrec_ZelditchMD.pdf")
  #Golden moles
  #pdf(file="mands_trc+gmole_gmole_ZelditchMD.pdf")
  
  #Non-microgale tenrecs
  #pdf(file="mands_nonmictrc+gmole_nonmic_tenrec_ZelditchMD.pdf")
#--------------------------------------------------------------------
par(mgp=c(3.2,0.5,0))
par(mar=c(5,7,4,2)+0.1) 

  md.hist <- hist(sim.md, xlab="ZelditchMD", main=NULL, las=1,
            xlim=c(obs.md, max(sim.md)), cex.lab=1.2)
    arrow.to.x.point(md.hist, obs.md, fraction.of.yaxis=50, line.fraction.of.yaxis=4,
                    height.above.xaxis=5, head.length=0.15, colour="blue", line.width=2.5)
  dev.off()

