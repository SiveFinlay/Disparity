#12/05/2014
#General script for simulating shape evolution across phylogenies and calculating disparity
    #Based on the previous script (30_04_14) but tidier because I'm using my new general functions

#Steps
  #1) Read in phylogenies, shape data and taxonomy
  #2) Choose which family (tenrecs or golden moles) 
  #3) Prune phylogenies
  #4) Shape simulation across phylogenies
  #5) PCA analysis of each simulation
  #6) Calculate disparity for each simulation and observed data
  #7) Compare observed and simulated disparity
  #8) Create output files
  
#########################################################
library(ape)
library(geiger)
library(geomorph)

#-------------------------------------------------------------------------------------
#First option: working directories

#Run on my computer
	source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
	source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PValueFunction_FromDistribution.r")
  source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r" )

#On the alien: save everything onto the USB
#  setwd("E:/Disparity")
#  source("E:/Disparity/DisparityFunctions_Variance_Range.r")
#  source("E:/Disparity/PValueFunction_FromDistribution.r")

######################################################
#1) READ IN DATA
######################################################

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
  
  #2b) Non-microgale tenrecs and all golden moles
      #shape coordinates
#      sps.mean <- dget(file="SkDors_nonmictenrec+gmole_sps.mean.txt")
      #taxonomic information
#      tax <- read.table("SkDors_nonmictenrec+gmole_taxonomy.txt")
#------------------------------------------------------
#SkLat

#1) Phylogenies
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
#     mytrees <- read.tree("SkLat_tenrec+gmole_101trees.phy")
     
#2) Data
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/sklat")
  #2a) All tenrecs and golden moles
      #shape coordinates
#      sps.mean <- dget(file="SkLat_tenrec+gmole_sps.mean.txt")
      #taxonomic information
#      tax <- read.table("SkLat_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
#    sps.mean <- dget(file="SkLat_nonmictenrec+gmole_sps.mean.txt")
    #taxonomic information
#    tax <- read.table("SkLat_nonmictenrec+gmole_sps.mean_taxonomy.txt")
#------------------------------------------------------
#SkVent
#1) Phylogenies
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
#    mytrees <- read.tree("SkVent_tenrec+gmole_101trees.phy")

#2) Data
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skvent")
  #2a) All tenrecs and golden moles
      #shape coordinates
#      sps.mean <- dget(file="SkVent_tenrec+gmole_sps.mean.txt")
      #taxonomic information
#     tax <- read.table("SkVent_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
    #sps.mean <- dget("SkVent_nonmictenrec+gmole_sps.mean.txt")
    #taxonomic information
    #tax <- read.table("SkVent_nonmictenrec+gmole_sps.mean_taxonomy.txt")
#------------------------------------------------------
#Mandibles
#1) Phylogenies
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/phylogenies")
#    mytrees <- read.tree("Mands_tenrec+gmole_101trees.phy")

#2) Data
#setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")
  #2a) All tenrecs and golden moles
      #shape coordinates
#      sps.mean <- dget(file="Mands_tenrec+gmole_sps.mean.txt")
      #taxonomic information
#     tax <- read.table("Mands_tenrec+gmole_sps.mean_taxonomy.txt")
  
  #2b) Non-microgale tenrecs and all golden moles
    #shape coordinates
#    sps.mean <- dget(file="Mands_nonmictenrec+gmole_sps.mean.txt")
    #taxonomic information
#    tax <- read.table("Mands_nonmictenrec+gmole_sps.mean_taxonomy.txt")

#################################################
#2) CHOOSE WHICH FAMILY 
#################################################
#Golden moles
# sps.tax <- tax$Binomial[which(tax$Family == "Chrysochloridae")]

#Tenrecs
  sps.tax <- tax$Binomial[which(tax$Family == "Tenrecidae")]

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

#Separate variance covariance matrix of the shape data for each of the phylogenies
  varcov <- as.list(rep(NA,length(sps.trees)))
    for(i in 1:length(sps.trees)){
      varcov[[i]] <- vcv(phy=sps.trees[[i]],twoDshape)
    }   

#simulate shape evolution on each phylogeny
  shape.sim <- as.list(rep(NA,length(sps.trees)))
    for (i in 1: length(sps.trees)){
      shape.sim[[i]] <- sim.char(sps.trees[[i]], varcov[[i]],nsim=1000,model="BM")
    }
    
#Combine simulations into one list
simlist <- list.matrix.to.array(shape.sim)  
######################################################
#5) PCA ANALYSIS
######################################################

#a) Simulated data
  shape.simPC <- calc.each.list(mylist=simlist, calculation=prcomp)

#Select the PC axes that account for 95% of the variation
  shape.simPC95axes <- NULL
    for (i in 1: length(shape.simPC)){
      #if PC1 axis explains less than 0.956% of the cumulative variation
      if(summary(shape.simPC[[i]])$importance[3,1] <= 0.956){
        #select the number of axes which do explain less than or equal to 0.956% of the total variation
        shape.simPC95axes[[i]] <- which(summary(shape.simPC[[i]])$importance[3,] <= 0.956)
      } else{
     #otherwise just select the first PC axis
      #NB: this could mean either that PC1 explains more than 95% of the variation
      #or that PC1 is less than 95% and PC2 is greater than 95%
        shape.simPC95axes[[i]] <- 1
        }
  }
 

#Use these numbers to select the corresponding PC axes for each simulation
  shape.simPC95 <- NULL
    for(i in 1: length(shape.simPC95axes)){
        shape.simPC95[[i]] <- shape.simPC[[i]]$x[,shape.simPC95axes[[i]]]
    }
  
#Turn all of the elements into matrices -makes it easier for later calculations
  shape.simPC95 <- calc.each.list(mylist=shape.simPC95, calculation=as.matrix) 

#b) Observed data
  #Do a principal components analysis on the family's (tenrec or golden mole)shape values only
    #i.e. don't use the PC axes from the global principal components analysis of all the species together
    #Makes sense because simulations are based on the shape coordinates for one family only

#PCA of mean shape values for each species
  obsPC <- prcomp(twoDshape)

#select the PC axes which correspond to 95% of the variation
  obsPC95 <- obsPC$x[,which(summary(obsPC)$importance[3,]<=0.956)]

##########################################################
#6) CALCULATE DISPARITY FOR SIMULATIONS AND OBSERVED DATA
##########################################################
#Simulated data
  #Variance measures
    sumvar <- calc.each.list(mylist=shape.simPC95, calculation=PCsumvar)
    prodvar <- calc.each.list(mylist=shape.simPC95, calculation=PCprodvar)
    
      #Matrix of the sum and product of variance for each simulation
      simPC.var <- matrix(NA,nrow=length(sumvar),ncol=2)
        colnames(simPC.var) <- c("SumVar","ProdVar")
        simPC.var[,1] <- unlist(sumvar)
        simPC.var[,2] <- unlist(prodvar)
  

  #Range measures
    sumrange <- calc.each.list(mylist=shape.simPC95, calculation=PCsumrange)
    prodrange <- calc.each.list(mylist=shape.simPC95, calculation=PCprodrange)
      
      #Matrix of the sum and product of range for each simulation
      simPC.range <- matrix(NA,nrow=length(sumrange),ncol=2)
        colnames(simPC.range) <- c("SumRange","ProdRange")
        simPC.range[,1] <- sumrange
        simPC.range[,2] <- prodrange
    
    
#Observed data
  #Variance
    obssumvar<-PCsumvar(obsPC95)
    obsprodvar<-PCprodvar(obsPC95)
  
  #Range
    obssumrange<-PCsumrange(obsPC95)
    obsprodrange<-PCprodrange(obsPC95)

######################################################### 
#7)COMPARE OBSERVED AND SIMULATED DISPARITY
#########################################################

#Compare observed disparity to the distribution of simulated values
  # (histograms in the output section below)
  sumvar.p<-pvalue.dist(distribution=simPC.var[,1], obs.val=obssumvar)
  prodvar.p<-pvalue.dist(distribution=simPC.var[,2], obs.val=obsprodvar)
  sumrange.p<-pvalue.dist(distribution=simPC.range[,1], obs.val=obssumrange)
  prodrange.p<-pvalue.dist(distribution=simPC.range[,2], obs.val=obsprodrange)



#Create atable to compare the disparity measures
  disp<-matrix(NA, nrow=4, ncol=5)
  rownames(disp)<-c("SumVar","ProdVar","SumRange","ProdRange")
  colnames(disp)<-c("Observed","Sim.min","Sim.max", "Sdev.sim","p.value")

    disp[1,1]<-obssumvar
    disp[1,2]<-range(simPC.var[,1])[1]
    disp[1,3]<-range(simPC.var[,1])[2]
    disp[1,4]<-sd(simPC.var[,1])
    disp[1,5]<-sumvar.p

    disp[2,1]<-obsprodvar
    disp[2,2]<-range(simPC.var[,2])[1]
    disp[2,3]<-range(simPC.var[,2])[2]
    disp[2,4]<-sd(simPC.var[,2])
    disp[2,5]<-prodvar.p

    disp[3,1]<-obssumrange
    disp[3,2]<-range(simPC.range[,1])[1]
    disp[3,3]<-range(simPC.range[,1])[2]
    disp[3,4]<-sd(simPC.range[,1])
    disp[3,5]<-sumrange.p

    disp[4,1]<-obsprodrange
    disp[4,2]<-range(simPC.range[,2])[1]
    disp[4,3]<-range(simPC.range[,2])[2]
    disp[4,4]<-sd(simPC.range[,2])
    disp[4,5]<-prodrange.p

#######################################
#8) CREATE THE OUTPUT FILES
#######################################
  #Change output folder and file names depending on the data

#a) Tables
  

#SkDors
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/output/shape_simulation/skdors")
  #Full data
#      write.table(file="SkDors_trc+gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkDors_trc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #Non-Microgale tenrecs
#    write.table(file="SkDors_nonmictrc+gmole_nonmictenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkDors_nonmictrc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#----------------------------------------------------------------------------------

#SkLat
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/output/shape_simulation/sklat")
  #Full data
#      write.table(file="SkLat_Trc+Gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkLat_Trc+Gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #Non-Microgale tenrecs
#      write.table(file="SkLat_nonmictrc+gmole_nonmictenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkLat_nonmictrc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#---------------------------------------------------------------------------------

#SkVent
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/output/shape_simulation/skvent")
  #Full data
#      write.table(file="SkVent_Trc+Gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkVent_Trc+Gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #Non-Microgale tenrecs
#      write.table(file="SkVent_nonmictrc+gmole_nonmictenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="SkVent_nonmictrc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#---------------------------------------------------------------------------------

#Mandibles
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/output/shape_simulation/mands")
  #Full data
#      write.table(file="Mands_Trc+Gmole_tenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="Mands_Trc+Gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #Non-Microgale tenrecs
#      write.table(file="Mands_nonmictrc+gmole_nonmictenrec_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#      write.table(file="Mands_nonmictrc+gmole_gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)

#**************************************
#b) PDF plots of the distribution of simulated disparity measures


#SkDors
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/SkDors_disparity_results")
#Working directory for the alien
#    setwd("E:/Disparity/SkDors_disparity_results")
  #Full data
#    pdf(file="SkDors_Trc+Gmole_TenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkDors_Trc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
  #Non-Microgale tenrecs
#    pdf(file="SkDors_NonMicTrc+Gmole_NonMicTenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkDors_NonMicTrc+Gmole_GmoleDisparity_1000simsDistribution.pdf")

#SkLat
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/SkLat_disparity_results")
#Working directory for the alien
#    setwd("E:/Disparity/SkLat_disparity_results")
  #Full data
#    pdf(file="SkLat_Trc+Gmole_TenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkLat_Trc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
  #Non-Microgale tenrecs
#    pdf(file="SkLat_NonMicTrc+Gmole_NonMicTenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkLat_NonMicTrc+Gmole_GmoleDisparity_1000simsDistribution.pdf")

#SkVent
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/SkVent_disparity_results")
#Working directory for the alien
#    setwd("E:/Disparity/SkVent_disparity_results")
  #Full data
#    pdf(file="SkVent_Trc+Gmole_TenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkVent_Trc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
  #Non-Microgale tenrecs
#    pdf(file="SkVent_NonMicTrc+Gmole_NonMicTenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="SkVent_NonMicTrc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
    
#Mandibles
#  setwd("C:/Users/sfinlay/Desktop/R_work/Morphometrics/Disparity/Mandibles_disparity_results")
#Working directory for the alien
#    setwd("E:/Disparity/Mandibles_disparity_results")
  #Full data
#    pdf(file="Mands_Trc+Gmole_TenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="Mands_Trc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
  #Non-Microgale tenrecs
#   pdf(file="Mands_NonMicTrc+Gmole_NonMicTenrecDisparity_1000simsDistribution.pdf")
#    pdf(file="Mands_NonMicTrc+Gmole_GmoleDisparity_1000simsDistribution.pdf")
    
par(mfrow=c(2,2))

hist(simPC.var[,1], xlab="SumVar",main=NULL)
hist(simPC.var[,2], xlab="SumProd",main=NULL)
hist(simPC.range[,1], xlab="SumRange",main=NULL)
hist(simPC.range[,2], xlab="ProdRange",main=NULL)

#Uncomment this line when I'm saving the plots above
#dev.off()


sv.hist<-hist(simPC.var[,1], xlab="SumVar", main=NULL)
obs.sv<-c(obssumvar,(mean(sv.hist$counts)/2))


hist(simPC.var[,1], xlab="SumVar", main=NULL)
points(obs.sv[1], obs.sv[2], pch=19, col="blue", cex=1.3)
arrows((obs.sv[1]),(obs.sv[2]+40),(obs.sv[1]),(obs.sv[2]+5), length=0.2, col="blue", lwd=2.5, code=2) 
#------------------------------------------
x<-rnorm(1000)
norm.hist<-hist(x)
p.hist<-c(-1.75, (mean(norm.hist$counts)/2))

hist(x)
points(p.hist[1], p.hist[2], pch=19, col="blue", cex=1.3)

#maybe change to density histogram plots?
#write functions for adding the arrows, distance of the arrow as a percentage of the y axis





