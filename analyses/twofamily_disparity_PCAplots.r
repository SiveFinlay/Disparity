#25/07/2014
  #Script to make a 4-panel PCA plot of each data set for my paper
  
  #Analysis code is exactly the same as the twofamily_disparity script except that the
  #data objects are re-named so I can have all four data sets within the one R session to plot the graphs together

#-----------------------------------------------------------------------------------

library(geomorph)
library(vegan)
library(boot)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")
########################################################
#READ IN DATA; directory will change for each data set
  #Re-name the objects with skdors, skvent etc to distinguish between them
########################################################
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data")
#SkDors data

#1) Landmarks
#landmarks + curves file with the control lines removed
  skdors.land <- readland.tps(file="skdors/Skdors_16_12_13_10landmarks+4curves_edited.TPS")

#2) Sliders
#edited sliders file (top 2 rows removed and the words before slide after put in instead
  skdors.curves <- as.matrix(read.table("skdors/Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
#file that has the correct taxonomy for each of the images
  skdors.taxa <- read.csv ("skdors/Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)

#4) Specimens to remove
  #Null
#--------------------------------------------------------
#SkLat data

#1) Landmarks
  sklat.land <- readland.tps(file="sklat/SkLat_08_11_13_9landmarks_2curves_edited.TPS")
#2) Sliders
  sklat.curves <- as.matrix(read.table(file="sklat/SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
#3) Taxonomy
  sklat.taxa <- read.csv("sklat/SkLat_08_11_13_Specimens+images.csv", header=TRUE)
#4) Specimens to remove
  sklat.rem <- read.csv("sklat/SkLat_remove_spec.csv", header=T)

#-----------------------------------------------------
#SkVent data

#1) Landmarks
  skvent.land <- readland.tps(file="skvent/SkVent_30_10_13_13landmarks+1curve_edited.TPS")
#2) Sliders
  skvent.curves <- as.matrix(read.table(file="skvent/SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
#3) Taxonomy
  skvent.taxa <- read.csv("skvent/SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
#4) Specimens to remove
  skvent.rem <- read.csv("skvent/SkVent_remove_spec.csv", header=T)
#------------------------------------------
#Mandibles data

#1) Landmarks
  mands.land <- readland.tps(file="mands/Mands_14_03_2014_7landmarks+4curves_edited.TPS")
#2) Sliders
  mands.curves <- as.matrix(read.table("mands/Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
#3) Taxonomy
  mands.taxa <- read.csv("mands/Mands_14_03_2014_Images+Specimens.csv", header=T)
#4) Specimens to remove
  mands.rem <- read.csv("mands/Mands_remove_spec.csv", header=T)

#################################################
#CLEAN UP THE DATA
#################################################
#Combine the taxonomic information with the arrays of landmark data

  #SkDors
    skdors.combine <- list(land=skdors.land, curves=skdors.curves, ID=skdors.taxa$ID,SpecID=skdors.taxa$SpecID, Order=skdors.taxa$Order_05,
                      Fam=skdors.taxa$Family_05, Genus=skdors.taxa$Genus_05, Species=skdors.taxa$Species_05, Binom=skdors.taxa$Binomial_05)

  # Remove the _sp. specimens
    skdors.sp <- which(skdors.combine$Species=="sp.")
    skdors.combine <- remove.from.list(skdors.combine, skdors.sp)
    skdors.combine <- droplevels.from.list(skdors.combine)
    
  # No specimens with missing data in skdors

  # Select the tenrec and golden mole specimens only
    skdors.tc.gm <- c(which(skdors.combine$Fam=="Chrysochloridae"), which(skdors.combine$Fam=="Tenrecidae"))
    skdors.data <- select.from.list(skdors.combine, skdors.tc.gm)
    skdors.data <- droplevels.from.list(skdors.data)
#----------------------------------------------
  #SkVent
    skvent.combine <- list(land=skvent.land, curves=skvent.curves, ID=skvent.taxa$ID, SpecID=skvent.taxa$SpecID, Order=skvent.taxa$Order_05,
                      Fam=skvent.taxa$Family_05, Genus=skvent.taxa$Genus_05, Species=skvent.taxa$Species_05, Binom=skvent.taxa$Binomial_05)

  # Remove the _sp. specimens
    skvent.sp <- which(skvent.combine$Species=="sp.")
    skvent.combine <- remove.from.list(skvent.combine, skvent.sp)
    skvent.combine <- droplevels.from.list(skvent.combine)

  # Remove specimens with missing data
    skvent.matching <- matching.id(skvent.rem$SpecID, skvent.combine$SpecID)
    skvent.combine <- remove.from.list(skvent.combine, skvent.matching)
    skvent.combine <- droplevels.from.list(skvent.combine)
    
  # Select the tenrec and golden mole specimens only
    skvent.tc.gm <- c(which(skvent.combine$Fam=="Chrysochloridae"), which(skvent.combine$Fam=="Tenrecidae"))
    skvent.data <- select.from.list(skvent.combine, skvent.tc.gm)
    skvent.data <- droplevels.from.list(skvent.data)
#----------------------------------------------
  #SkLat
    sklat.combine <- list(land=sklat.land, curves=sklat.curves, ID=sklat.taxa$ID, SpecID=sklat.taxa$SpecID, Order=sklat.taxa$Order_05,
                      Fam=sklat.taxa$Family_05, Genus=sklat.taxa$Genus_05, Species=sklat.taxa$Species_05, Binom=sklat.taxa$Binomial_05)

  # Remove the _sp. specimens
    sklat.sp <- which(sklat.combine$Species=="sp.")
    sklat.combine <- remove.from.list(sklat.combine, sklat.sp)
    sklat.combine <- droplevels.from.list(sklat.combine)

  # Remove specimens with missing data
    sklat.matching <- matching.id(sklat.rem$SpecID, sklat.combine$SpecID)
    sklat.combine <- remove.from.list(sklat.combine, sklat.matching)
    sklat.combine <- droplevels.from.list(sklat.combine)
    
  # Select the tenrec and golden mole specimens only
    sklat.tc.gm <- c(which(sklat.combine$Fam=="Chrysochloridae"), which(sklat.combine$Fam=="Tenrecidae"))
    sklat.data <- select.from.list(sklat.combine, sklat.tc.gm)
    sklat.data <- droplevels.from.list(sklat.data)
#----------------------------------------------
  #Mandibles
    mands.combine <- list(land=mands.land, curves=mands.curves, ID=mands.taxa$ID, SpecID=mands.taxa$SpecID, Order=mands.taxa$Order_05,
                      Fam=mands.taxa$Family_05, Genus=mands.taxa$Genus_05, Species=mands.taxa$Species_05, Binom=mands.taxa$Binomial_05)

  # Remove the _sp. specimens
    mands.sp <- which(mands.combine$Species=="sp.")
    mands.combine <- remove.from.list(mands.combine, sklat.sp)
    mands.combine <- droplevels.from.list(mands.combine)

  # Remove specimens with missing data
    mands.matching <- matching.id(mands.rem$SpecID, mands.combine$SpecID)
    mands.combine <- remove.from.list(mands.combine, mands.matching)
    mands.combine <- droplevels.from.list(mands.combine)

  # Select the tenrec and golden mole specimens only
    mands.tc.gm <- c(which(mands.combine$Fam=="Chrysochloridae"), which(mands.combine$Fam=="Tenrecidae"))
    mands.data <- select.from.list(mands.combine, mands.tc.gm)
    mands.data <- droplevels.from.list(mands.data)


#######################################
#PROCRUSTES SUPERIMPOSTION
#######################################

#SkDors
  #General Procrustes Alignment
    skdors.dataGPA <- gpagen(skdors.data$land, curves=skdors.data$curves, ProcD=TRUE,)

  #List the coordinates with the taxonomic information
    skdors.proc <- list(coords=skdors.dataGPA$coords, csize=skdors.dataGPA$Csize, ID=skdors.data$ID, SpecID=skdors.data$SpecID,
                        Order=skdors.data$Order, Fam=skdors.data$Fam, Genus=skdors.data$Genus, Species=skdors.data$Species, Binom=skdors.data$Binom)

#SkVent
    skvent.dataGPA <- gpagen(skvent.data$land, curves=skvent.data$curves, ProcD=TRUE,)

    skvent.proc <- list(coords=skvent.dataGPA$coords, csize=skvent.dataGPA$Csize, ID=skvent.data$ID, SpecID=skvent.data$SpecID,
                        Order=skvent.data$Order, Fam=skvent.data$Fam, Genus=skvent.data$Genus, Species=skvent.data$Species, Binom=skvent.data$Binom)

#SkLat
    sklat.dataGPA <- gpagen(sklat.data$land, curves=sklat.data$curves, ProcD=TRUE,)

    sklat.proc <- list(coords=sklat.dataGPA$coords, csize=sklat.dataGPA$Csize, ID=sklat.data$ID, SpecID=sklat.data$SpecID,
                        Order=sklat.data$Order, Fam=sklat.data$Fam, Genus=sklat.data$Genus, Species=sklat.data$Species, Binom=sklat.data$Binom)

#Mands
    mands.dataGPA <- gpagen(mands.data$land, curves=mands.data$curves, ProcD=TRUE,)
    
    mands.proc <- list(coords=mands.dataGPA$coords, csize=mands.dataGPA$Csize, ID=mands.data$ID, SpecID=mands.data$SpecID,
                        Order=mands.data$Order, Fam=mands.data$Fam, Genus=mands.data$Genus, Species=mands.data$Species, Binom=mands.data$Binom)

#######################################
#SPECIES AVERAGING
#######################################
#SkDors
  #Group the arrays of coordinates according to species
    skdors.group.sps.coords <- species.coordinates(skdors.proc$coords, skdors.proc$Binom)

  #Average coordinate values for each species
    skdors.sps.mean <- mean.coords(skdors.group.sps.coords)

#SkVent
    skvent.group.sps.coords <- species.coordinates(skvent.proc$coords, skvent.proc$Binom)
    skvent.sps.mean <- mean.coords(skvent.group.sps.coords)

#SkLat
    sklat.group.sps.coords <- species.coordinates(sklat.proc$coords, sklat.proc$Binom)
    sklat.sps.mean <- mean.coords(sklat.group.sps.coords)

#Mands
    mands.group.sps.coords <- species.coordinates(mands.proc$coords, mands.proc$Binom)
    mands.sps.mean <- mean.coords(mands.group.sps.coords)


#######################################
#PRINCIPAL COMPONENTS ANALYSIS
#######################################

#SkDors
  skdors.sps.meanPCA <- plotTangentSpace(skdors.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

  #First and second PC axes
    skdors.xaxis <- skdors.sps.meanPCA$x[,1]
    skdors.yaxis <- skdors.sps.meanPCA$x[,2]

  #Data frame of the unique family and binomical combinations; use it to colour points by family
    skdors.sp.fam <- as.data.frame(unique(cbind(as.matrix(skdors.proc$Fam), as.matrix(skdors.proc$Binom))))
    colnames(skdors.sp.fam) <- c("Family","Binomial")
#-------------------------------------
#SkVent
  skvent.sps.meanPCA <- plotTangentSpace(skvent.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

  #First and second PC axes
    skvent.xaxis <- skvent.sps.meanPCA$x[,1]
    skvent.yaxis <- skvent.sps.meanPCA$x[,2]

  #Data frame of the unique family and binomical combinations; use it to colour points by family
    skvent.sp.fam <- as.data.frame(unique(cbind(as.matrix(skvent.proc$Fam), as.matrix(skvent.proc$Binom))))
    colnames(skvent.sp.fam) <- c("Family","Binomial")
#------------------------------------------
#SkLat
  sklat.sps.meanPCA <- plotTangentSpace(sklat.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

  #First and second PC axes
    sklat.xaxis <- sklat.sps.meanPCA$x[,1]
    sklat.yaxis <- sklat.sps.meanPCA$x[,2]

  #Data frame of the unique family and binomical combinations; use it to colour points by family
    sklat.sp.fam <- as.data.frame(unique(cbind(as.matrix(sklat.proc$Fam), as.matrix(sklat.proc$Binom))))
    colnames(sklat.sp.fam) <- c("Family","Binomial")
#-------------------------------------------
#Mandibles
  mands.sps.meanPCA <- plotTangentSpace(mands.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

  #First and second PC axes
    mands.xaxis <- mands.sps.meanPCA$x[,1]
    mands.yaxis <- mands.sps.meanPCA$x[,2]

  #Data frame of the unique family and binomical combinations; use it to colour points by family
    mands.sp.fam <- as.data.frame(unique(cbind(as.matrix(mands.proc$Fam), as.matrix(mands.proc$Binom))))
    colnames(mands.sp.fam) <- c("Family","Binomial")

#--------------------------------------
#PCA plots: default colour palette so Chrysochloridae are black and Tenrecidae are red

par(mfrow=c(2,2))
#SkDors
  plot(skdors.xaxis, skdors.yaxis, xlab="", ylab="", las=1,
       col=skdors.sp.fam$Family, pch=16, bty="l", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1, at=c(round(min(skdors.xaxis),3), 0, round(max(skdors.xaxis),3)), las=1, cex=1.5)
    #same for the y axis
      axis(side=2, at=c(round(min(skdors.yaxis),3), 0, round(max(skdors.yaxis),3)), las=1, cex=1.5)
    #add dotted lines along 0,0
      abline(0,0, h=0, v=0, lty=2, lwd=1)
      
#SkVent
  plot(skvent.xaxis, skvent.yaxis, xlab="", ylab="", las=1,
       col=skvent.sp.fam$Family, pch=16, bty="l", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1, at=c(round(min(skvent.xaxis),3), 0, round(max(skvent.xaxis),3)), las=1, cex=1.5)
    #same for the y axis
      axis(side=2, at=c(round(min(skvent.yaxis),3), 0, round(max(skvent.yaxis),3)), las=1, cex=1.5)
    #add dotted lines along 0,0
      abline(0,0, h=0, v=0, lty=2, lwd=1)


#SkLat
  plot(sklat.xaxis, sklat.yaxis, xlab="", ylab="", las=1,
       col=sklat.sp.fam$Family, pch=16, bty="l", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1, at=c(round(min(sklat.xaxis),3), 0, round(max(sklat.xaxis),3)), las=1, cex=1.5)
    #same for the y axis
      axis(side=2, at=c(round(min(sklat.yaxis),3), 0, round(max(sklat.yaxis),3)), las=1, cex=1.5)
    #add dotted lines along 0,0
      abline(0,0, h=0, v=0, lty=2, lwd=1)
      

#Mands
  plot(mands.xaxis, mands.yaxis, xlab="", ylab="", las=1,
       col=mands.sp.fam$Family, pch=16, bty="l", cex.lab=1, cex=1.2, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1, at=c(round(min(mands.xaxis),3), 0, round(max(mands.xaxis),3)), las=1, cex=1.5)
    #same for the y axis
      axis(side=2, at=c(round(min(mands.yaxis),3), 0, round(max(mands.yaxis),3)), las=1, cex=1.5)
    #add dotted lines along 0,0
      abline(0,0, h=0, v=0, lty=2, lwd=1)

#I copied and pasted the plots into a powerpoint presentation and added the labels to the axes

#--------------------------------------
#Select PC axes for each analysis
  #I need these numbers as a comparison of the number of axes that I used for each data set
  skdors.PC95axes <- selectPCaxes(skdors.sps.meanPCA, 0.956, skdors.sps.mean$Binom)
  skvent.PC95axes <- selectPCaxes(skvent.sps.meanPCA, 0.956, skvent.sps.mean$Binom)
  sklat.PC95axes <- selectPCaxes(sklat.sps.meanPCA, 0.956, sklat.sps.mean$Binom)
  mands.PC95axes <- selectPCaxes(mands.sps.meanPCA, 0.956, mands.sps.mean$Binom)