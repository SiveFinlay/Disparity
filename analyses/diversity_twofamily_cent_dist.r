#16/09/14

  #Simplified disparity calculations
    #Calculate disparity as based on euclidean distances from the centroid for each family
      #Compare the dispersion of tenrecs and golden moles around their average skull shapes

  #Overall result for skulls: no significant difference for the full data
    #Tenrecs>golden moles in the reduced analysis with the Microgale subsample
    
    #Leaving out the mandibles because they have G>T for all data and no significant difference with the subsample data set

#steps:
  #1) Read in a clean up raw landmark data
    #OPTIONS; choices depending on the analysis
  #2) Procrustes superimposition of tenrecs and golden moles
  #3) Find the average Procrustes shape coordinates for each species
  #4) PCA of the shape coordinates
  #5) Select PC axes that account for 95% of the variation
  #6) Calculate disparity based on distances from the centroid

library(geomorph)
library(vegan)
library(plotrix)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Morpho_diversity_functions.r")

setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/")

#######################
#Read in the data: 3 views of skulls, 2 options for mandibles (1 or 4 curves)
#################
#SkDors
    #1) Landmarks
      #land <- readland.tps(file="skdors/Skdors_16_12_13_10landmarks+4curves_edited.TPS")
    #2) Sliders
      #curves <- as.matrix(read.table("skdors/Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))
    #3) Taxonomy
      #taxa <- read.csv ("skdors/Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)
    #4) Specimens to remove
      #Null

#-----------------------------------------------------
#SkVent
  #1) Landmarks
    #land <- readland.tps(file="skvent/SkVent_30_10_13_13landmarks+1curve_edited.TPS")
  #2) Sliders
    #curves <- as.matrix(read.table(file="skvent/SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
  #3) Taxonomy
    #taxa <- read.csv("skvent/SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
  #4) Specimens to remove
    #rem <- read.csv("skvent/SkVent_remove_spec.csv", header=T)

#--------------------------------------------------------
#SkLat
  #1) Landmarks
    land <- readland.tps(file="sklat/SkLat_08_11_13_9landmarks_2curves_edited.TPS")
  #2) Sliders
    curves <- as.matrix(read.table(file="sklat/SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
  #3) Taxonomy
    taxa <- read.csv("sklat/SkLat_08_11_13_Specimens+images.csv", header=TRUE)
  #4) Specimens to remove
    rem <- read.csv("sklat/SkLat_remove_spec.csv", header=T)



#################################################
#CLEAN UP THE DATA
#################################################
#Combine the taxonomic information with the array of landmark data
  combine <- list(land=land, curves=curves, ID=taxa$ID,SpecID=taxa$SpecID, Order=taxa$Order_05,
                  Fam=taxa$Family_05, Genus=taxa$Genus_05, Species=taxa$Species_05, Binom=taxa$Binomial_05)

# Remove the _sp. specimens
 sp <- which(combine$Species=="sp.")

 combine <- remove.from.list(combine, sp)
 combine <- droplevels.from.list(combine)

#********************************************
#OPTION; depending on the data and the analysis
#************************************
#Remove the specimens listed in rem (sklat, skvent and mands data)
  #doesn't apply to the skdors data because rem is NULL

#find the ID numbers of specimens with missing data
  matching <- matching.id(rem$SpecID, combine$SpecID)
    combine <- remove.from.list(combine, matching)
    combine <- droplevels.from.list(combine)

##################################
#Select specimens that you want
#####################################
#Select the tenrecs and golden moles only
  tc.gm <- c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))

  mydata <- select.from.list(combine, tc.gm)
  mydata <- droplevels.from.list(mydata)

#************************************
#Option to select a subset of tenrecs
#************************************
#I originally removed all of the Microgale but it makes more sense to keep at least some of them

#Find all of the rows that are Microgale specimens
   mic <- which(mydata$Genus=="Microgale")
#Find how many different Microgale species there are
  mic.data <- select.from.list(mydata, mic)
  mic.data <- droplevels.from.list(mic.data)

  #Soarimalala et al 2011 divide Microgale into 5 groups based on body size and tail length
    #I'm using these as proxies for diversity across the Microgale genus
      #Select 1 species to represent each of the 5 groups: parvula, brevicaudata, dryas, longicaudata, dobsoni

  #Row numbers for the selected microgale species
    sel.mic.id <- sort(c(which(mic.data$Species == "parvula"), which(mic.data$Species == "brevicaudata"),
                  which(mic.data$Species == "dryas"), which(mic.data$Species == "longicaudata"),
                  which(mic.data$Species == "dobsoni")))
                  
  #List of Microgale species which are not the selected ones
   mic.spec.rem <- droplevels((remove.from.list(mic.data, sel.mic.id))$Binom)

  #Remove these Microgale (12 species that are not the 5 selected ones)

  #Find the ID numbers of those species within the main data set
  mic.rem.id <- NULL
    for (i in 1:length(levels(mic.spec.rem))){
     mic.rem.id[[i]] <- which(mydata$Binom == levels(mic.spec.rem)[i])
   }

   #Remove those IDs from the data
    mydata <- remove.from.list(mydata, unlist(mic.rem.id))
    mydata <- droplevels.from.list(mydata)

#End of the option to remove some tenrecs
#*****************


#######################################
#PROCRUSTES SUPERIMPOSTION
#######################################

#General Procrustes Alignment of all of the scaled coordinates
  mydataGPA <- gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE,)
  #ProcD=TRUE means that the coordinates are aligned by procrustes distance rather than bending energy
      # which means that RWA is equivalent to PCA (Zelditch 2012, page 150)

#List the coordinates with the taxonomic information
  Proc.co <- list(coords=mydataGPA$coords,csize=mydataGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
                  Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

#######################################
#SPECIES AVERAGING
#######################################

#group the arrays of coordinates according to species
  group.sps.coords <- species.coordinates(Proc.co$coords, Proc.co$Binom)

#average coordinate values for each species
  sps.mean <- mean.coords(group.sps.coords)

#list of species
  binom <- sps.mean$Binom

#######################################
#PRINCIPAL COMPONENTS ANALYSIS
#######################################

sps.meanPCA <- plotTangentSpace(sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

#Re-create the PCA plot
  xaxis <- sps.meanPCA$x[,1]
  yaxis <- sps.meanPCA$x[,2]

#colour points by family
  #data frame of the unique family and binomical combinations
  sp.fam <- as.data.frame(unique(cbind(as.matrix(Proc.co$Fam), as.matrix(Proc.co$Binom))))
    colnames(sp.fam) <- c("Family","Binomial")

#I haven't made a pretty PCA graph from this script yet

#######################################
#SELECT PC AXES
#######################################
PC95axes <- selectPCaxes(sps.meanPCA, 0.956, binom)
#NB: results could change depending on the threshold set for the number of axes to use
#But the dimensionality of the two families is the same

#select the rows on those axes that correspond to each family
  gmolePC <- PC95axes[which(sp.fam$Family=="Chrysochloridae"),]
  tenrecPC <- PC95axes[which(sp.fam$Family=="Tenrecidae"),]
##########################
#Compare positions in morphospace
########################

#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
  PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")

#################
#Diversity of families based on distances to centroid
#################

#Distance from each species to that family's centroid
  gmole.cent.dist <- euc.dist.cent (gmolePC)
  tenrec.cent.dist <- euc.dist.cent (tenrecPC)
  
#Compare the distances to centroids in tenrecs and golden moles
  cent.dist <- matrix(nrow=nrow(PC95axes), ncol=2)
    colnames(cent.dist) <- c("dist", "group")
    cent.dist[,1] <- c(gmole.cent.dist, tenrec.cent.dist)
    cent.dist[,2] <- c(rep("gmole", length(gmole.cent.dist)), rep("tenrec", length(tenrec.cent.dist)))


#Compare the two groups with a t test
  comp.cent <- t.test(as.numeric(cent.dist[,1]) ~ cent.dist[,2])

  
#Mean and standard error of those distances from the centroid
  mean.se <- matrix(nrow=2, ncol=2)
    colnames(mean.se) <- c("mean", "se")
    rownames(mean.se) <- c("gmole", "tenrec")
  mean.se[,1] <- c(mean(gmole.cent.dist), mean(tenrec.cent.dist))
  mean.se[,2] <- c(std.error(gmole.cent.dist), std.error(tenrec.cent.dist))

################################
#Pairwise permutation tests
##############################
#Check that the significant differences are not just an artefact of differences in sample size

#Observed differences in the mean
  obs.mean.diff <- mean.se[2,1] - mean.se[1,1]  

#Permutation test for significant difference in the mean
  perm.mean <- group.diff(1000, sp.fam$Family, PC95axes, mean)
#Test for significant difference
  perm.mean.pvalue <- pvalue.dist(perm.mean, obs.mean.diff)


#Summary table of the results
      perm.res.summary <- matrix(NA, nrow=1, ncol=6)
        colnames(perm.res.summary) <- c("obs.tenrec", "obs.gmole", "obs.diff", "perm.min", "perm.max", "pvalue") 
        perm.res.summary[1,] <- c(mean.se[2,1], mean.se[1,1], obs.mean.diff, min(perm.mean), max(perm.mean), perm.mean.pvalue) 

#############################################
#PCA plot
#Families have different symbols instead of different colours

  #plot(xaxis, yaxis, xlab="", ylab="", las=1,
       #pch=c(16,17)[sp.fam$Family], bty="l", cex.lab=1.75, cex=1.5, xaxt="n", yaxt="n")
    #draw the min,max and 0 values on the x axis
      #axis(side=1, at=c(round(min(xaxis),3), 0, round(max(xaxis),3)), las=1, cex=1.75)
    #same for the y axis
      #axis(side=2, at=c(round(min(yaxis),3), 0, round(max(yaxis),3)), las=1, cex=1.75)
    #add dotted lines along 0,0
      #abline(0,0, h=0, v=0, lty=2, lwd=1)
##############################################
#Output

setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/")      
      
#Skdors: full data
  # PC.man
    #capture.output(PC.man, file= "skdors/morpho_diversity/skdors_manova.res.txt")
  # mean.se
    #capture.output(mean.se, file= "skdors/morpho_diversity/skdors_mean+se.txt")    
  # t test results
    #capture.output(comp.cent, file= "skdors/morpho_diversity/skdors_ttest.txt") 
  #perm.res.summary
    #capture.output(perm.res.summary, file= "skdors/morpho_diversity/skdors_perm.res.summary.txt") 
      
#Skdors: subset of tenrecs
  # PC.man
    #capture.output(PC.man, file= "skdors/morpho_diversity/skdors_subtenrec_manova.res.txt")
  # mean.se
    #capture.output(mean.se, file= "skdors/morpho_diversity/skdors_subtenrec_mean+se.txt")    
  # t test results
    #capture.output(comp.cent, file= "skdors/morpho_diversity/skdors_subtenrec_ttest.txt") 
  #perm.res.summary
    #capture.output(perm.res.summary, file= "skdors/morpho_diversity/skdors_subtenrec_perm.res.summary.txt") 
#--------------------------------------------------------------------------------------------
#Skvent: full data
  # PC.man
    #capture.output(PC.man, file= "skvent/morpho_diversity/skvent_manova.res.txt")
  # mean.se
    #capture.output(mean.se, file= "skvent/morpho_diversity/skvent_mean+se.txt")    
  # t test results
    #capture.output(comp.cent, file= "skvent/morpho_diversity/skvent_ttest.txt") 
  #perm.res.summary
    #capture.output(perm.res.summary, file= "skvent/morpho_diversity/skvent_perm.res.summary.txt") 
      
#Skvent: subset of tenrecs
  # PC.man
    #capture.output(PC.man, file= "skvent/morpho_diversity/skvent_subtenrec_manova.res.txt")
  # mean.se
    #capture.output(mean.se, file= "skvent/morpho_diversity/skvent_subtenrec_mean+se.txt")    
  # t test results
    #capture.output(comp.cent, file= "skvent/morpho_diversity/skvent_subtenrec_ttest.txt") 
  #perm.res.summary
    #capture.output(perm.res.summary, file= "skvent/morpho_diversity/skvent_subtenrec_perm.res.summary.txt") 
#--------------------------------------------------------------------------------------------
#Sklat: full data
  # PC.man
    #capture.output(PC.man, file= "sklat/morpho_diversity/sklat_manova.res.txt")
  # mean.se
    #capture.output(mean.se, file= "sklat/morpho_diversity/sklat_mean+se.txt")    
  # t test results
    #capture.output(comp.cent, file= "sklat/morpho_diversity/sklat_ttest.txt") 
  #perm.res.summary
    #capture.output(perm.res.summary, file= "sklat/morpho_diversity/sklat_perm.res.summary.txt") 
      
#Sklat: subset of tenrecs
  # PC.man
    capture.output(PC.man, file= "sklat/morpho_diversity/sklat_subtenrec_manova.res.txt")
  # mean.se
    capture.output(mean.se, file= "sklat/morpho_diversity/sklat_subtenrec_mean+se.txt")    
  # t test results
    capture.output(comp.cent, file= "sklat/morpho_diversity/sklat_subtenrec_ttest.txt") 
  #perm.res.summary
    capture.output(perm.res.summary, file= "sklat/morpho_diversity/sklat_subtenrec_perm.res.summary.txt") 
###########################################################
    

#############################################