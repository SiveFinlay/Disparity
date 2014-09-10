#10/09/2014
#Revised disparity calculations after the steering committee meeting
#It's exactly the same data and analysis  approach as the twofamily_disparity script but the code is tidier
#The only difference is that I don't remove all of the Microgale specimens: keep a subset of 5 species

#This script can be used for two options
  #1) Compare disparity in all tenrecs to all golden moles
  #2) Compare dispartiy in a subsample of tenrecs (5 of the 17 Microgale) to all golden moles
#So far I've only used this script for option 2 (option 1 was covered in the previous script)

#steps:
  #1) Read in a clean up raw landmark data
    #OPTIONS; choices depending on the analysis
  #2) Procrustes superimposition of tenrecs and golden moles
  #3) Find the average Procrustes shape coordinates for each species
  #4) PCA of the shape coordinates
  #5) Select PC axes that account for 95% of the variation
  #6) Calculate disparity measures
  #7) Compare disparity in families; npMANOVA
  #8) Compare disparity measures directly: permutation tests
  #9) Output files: shape data, taxonomy, disparity, disparity comparisons


#output from this script
  #data
    #1) The average shape coordinates for each species of the GPA-aligned specimens (sps.mean)
    #2) The taxonomic information (Family and Binomial) for these shape coordinates (sp.fam)
        # (the binomial names are in the same order in each object)
    #3) Table of npMANOVA reults; based on distance matrix and PC axes              (manova.res)
    #4) Summary table of results from permutation tests for significant differences in disparity (perm.res.summary)

  #figures
    #1) PCA plots  - I haven't added the code to this script


library(geomorph)
library(vegan)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")

setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/")

#######################
#Read in the data
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
#------------------------------------------
#Mandibles
  #1) Landmarks
    #land <- readland.tps(file="mands/Mands_14_03_2014_7landmarks+4curves_edited.TPS")
  #2) Sliders
    #curves <- as.matrix(read.table("mands/Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
  #3) Taxonomy
    #taxa <- read.csv("mands/Mands_14_03_2014_Images+Specimens.csv", header=T)
  #4) Specimens to remove
    #rem <- read.csv("mands/Mands_remove_spec.csv", header=T)
    
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

#select the axes for each family
  gmolePC <- PC95axes[which(sp.fam$Family=="Chrysochloridae"),]
  tenrecPC <- PC95axes[which(sp.fam$Family=="Tenrecidae"),]
  
#######################################
#CALCULATE DISPARITY
#######################################
#Based on PC axes (original script had an additional ZelditchMD metric)
  #Tenrecs
    #variance
    tenrec.v <- PCvariance(tenrecPC)
    tenrec.sv <- PCsumvar(tenrecPC)
    tenrec.pv <- PCprodvar(tenrecPC)

    #range
    tenrec.r <- PCrange(tenrecPC)
    tenrec.sr <- PCsumrange(tenrecPC)
    tenrec.pr <- PCprodrange(tenrecPC)

  #Golden moles
    #variance
    gmole.v <- PCvariance(gmolePC)
    gmole.sv <- PCsumvar(gmolePC)
    gmole.pv <- PCprodvar(gmolePC)

    #range
    gmole.r <- PCrange(gmolePC)
    gmole.sr <- PCsumrange(gmolePC)
    gmole.pr <- PCprodrange(gmolePC)

#Put the disparity calculations into a single table
  disp <- matrix(NA,nrow=2, ncol=4)
    rownames(disp) <- c("Tenrec","Gmole")
    colnames(disp) <- c("SumVar","ProdVar","SumRange","ProdRange")
    disp[,1] <- c(tenrec.sv, gmole.sv)
    disp[,2] <- c(tenrec.pv, gmole.pv)
    disp[,3] <- c(tenrec.sr, gmole.sr)
    disp[,4] <- c(tenrec.pr, gmole.pr)

#######################################
#COMPARE FAMILIES
#######################################
#Compare morphospace occupation (not comparing disparity metrics directly)
        #Not very meaningful because it just shows that overall families occupy significantly different places in morphospace
          #i.e. it doesn't show the specific differences between particular families

#1) Distance matrix

#Euclidean distance matrix of all of the species
  Euc.dist <- as.matrix(dist(PC95axes, method="euclidean", diag=FALSE,upper=FALSE))

#npMANOVA of the distance matrix separated by family
  dist.man <- adonis(Euc.dist~sp.fam$Family, data=sp.fam, permutations=9999,method="euclidean")
    #extract the f, r2 and p values
      dist.man.frp <- anova.frp(dist.man)

#2) PC axes
#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
  PC.man <- adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")
    #extract the f, r2 and p values
      PC.man.frp <- anova.frp(PC.man)

#Compare the distance and PC npMANOVA results
  manova.res <- rbind(dist.man.frp, PC.man.frp)
  rownames(manova.res) <-c ("dist.man", "PC.man")
  
#-------------------------------------
#Test for significant differences in disparity (modified code from Steve Wang, email on 10/06/2014)
  #Advantage of this method is that it takes differences in sample size into account

  #Pairwise observed differences in disparity among the two family groups
    #Tenrec disparity - Gmole disparity
    
    obs.diff.sv <- disp[2,1] - disp[1,1]
    obs.diff.pv <- disp[2,2] - disp[1,2]
    obs.diff.sr <- disp[2,3] - disp[1,3]
    obs.diff.pr <- disp[2,4] - disp[1,4]

#Permutation tests for significant differences in disparity between tenrecs and golden moles
      tc.gm.perm.sv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCsumvar)
      tc.gm.perm.pv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCprodvar)
      tc.gm.perm.sr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCsumrange)
      tc.gm.perm.pr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCprodrange)

      #test for significant differences
          #(NB: Another way to look at the distribution is with a table showing when the permutated values are greater than the observed)
          #tc.gm.obs.dist.sv <- table(tc.gm.perm.sv >= obs.diff.sv[1,3])
      tc.gm.pvalue.sv <- pvalue.dist(tc.gm.perm.sv, obs.diff.sv)
      tc.gm.pvalue.pv <- pvalue.dist(tc.gm.perm.pv, obs.diff.pv)
      tc.gm.pvalue.sr <- pvalue.dist(tc.gm.perm.sr, obs.diff.sr)
      tc.gm.pvalue.pr <- pvalue.dist(tc.gm.perm.pr, obs.diff.pr)


#Summary table of the results
      perm.res.summary <- matrix(NA, nrow=4, ncol=7)
        colnames(perm.res.summary) <- c("metric", "obs.tenrec", "obs.gmole", "obs.diff", "perm.min", "perm.max", "pvalue")
        perm.res.summary[,1] <- c("sumvar", "prodvar", "sumrange", "prodrange")
          #tenrec disparity, rounded to 3 significant figures
        perm.res.summary[,2] <- signif(disp[1,],3)
           #gmole disparity, rounded to 3 significant figures
        perm.res.summary[,3] <- signif(disp[2,],3)
         #observed differences, everything rounded to 3 significant digits
        perm.res.summary[,4] <- c(signif(obs.diff.sv,3), signif(obs.diff.pv,3), signif(obs.diff.sr,3), signif(obs.diff.pr,3))
         #minimum permutated values, everything rounded to 3 significant digits
        perm.res.summary[,5] <- c(signif(min(tc.gm.perm.sv),3), signif(min(tc.gm.perm.pv),3), signif(min(tc.gm.perm.sr),3), signif(min(tc.gm.perm.pr),3))
          #maximum permutated values, everything rounded to 3 significant digits
        perm.res.summary[,6] <- c(signif(max(tc.gm.perm.sv),3), signif(max(tc.gm.perm.pv),3), signif(max(tc.gm.perm.sr),3), signif(max(tc.gm.perm.pr),3))
          #p values from comparing the observed differences to the distribution of permutated values
        perm.res.summary[,7] <- c(tc.gm.pvalue.sv, tc.gm.pvalue.pv, tc.gm.pvalue.sr, tc.gm.pvalue.pr)

  perm.res.summary <- as.data.frame(perm.res.summary)

#######################################
#OUTPUT FILES
#######################################

setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/")

#**************************************************
#1) Shape data, taxonomy , disparity measures, disparity comparisons
#****************************************************
#SkDors
#Subset of tenrecs (17 species including 5 Microgale) compared to golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="skdors/SkDors_submic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="skdors/SkDors_submic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="skdors/SkDors_submic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="skdors/SkDors_submic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#----------------------------------------
#SkLat
#Subset of tenrecs (17 species including 5 Microgale) compared to golden moles
  #1) Average shape coordinates
     dput(sps.mean, file="sklat/SkLat_submic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
     write.table(file="sklat/SkLat_submic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
     write.table(file="sklat/SkLat_submic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
     write.table(file="sklat/SkLat_submic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#--------------------------------------
#SkVent
#Subset of tenrecs (17 species including 5 Microgale) compared to golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="skvent/SkVent_submic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="skvent/SkVent_submic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="skvent/SkVent_submic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="skvent/SkVent_submic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#-----------------------------------------------
#Mands
#Subset of tenrecs (17 species including 5 Microgale) compared to golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="mands/Mands_submic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="mands/Mands_submic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="mands/Mands_submic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="mands/Mands_submic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
