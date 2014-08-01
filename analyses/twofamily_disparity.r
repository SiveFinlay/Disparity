#09/05/2014
  #updated on 11/06/2014; permutation tests for significant differences in disparity (instead of original rarefaction)
  #updated on 16/07/2014: tests of additional ideas from Uri at BES Nottingham
    #1) compare disparity of the other family groups; on the basis that golden moles are all fossorial so you 
      #would expect them to have lower morphological variation than the more ecologically diverse tenrecs
    #2) find a way to look at intraspecific variation

#Re-coding and re-analysis of my disparity project

#general script which can be used with any of the data sets (skdors, skvent, sklat, mands)
  #change the input files depending on which data I'm using
  
#compare disparity in tenrecs and golden moles only  # change this to comparing all of the families
  #additional option of selecting just microgale tenrecs

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
    #1) PCA plots

    
#-----------------------------------------------------------------------------------

library(geomorph)
library(vegan)
library(boot)

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")
########################################################
#READ IN DATA; directory will change for each data set
########################################################
#SkDors data
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/skdors")

#1) Landmarks
#landmarks + curves file with the control lines removed
  land <- readland.tps(file="Skdors_16_12_13_10landmarks+4curves_edited.TPS")

#2) Sliders
#edited sliders file (top 2 rows removed and the words before slide after put in instead
  curves <- as.matrix(read.table("Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
#file that has the correct taxonomy for each of the images
  taxa <- read.csv ("Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)

#4) Specimens to remove
  #Null
#--------------------------------------------------------
#SkLat data
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/sklat")

#1) Landmarks
  #land <- readland.tps(file="SkLat_08_11_13_9landmarks_2curves_edited.TPS")
#2) Sliders
  #curves <- as.matrix(read.table(file="SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
#3) Taxonomy
  #taxa <- read.csv("SkLat_08_11_13_Specimens+images.csv", header=TRUE)
#4) Specimens to remove
  #rem <- read.csv("SkLat_remove_spec.csv", header=T)

#-----------------------------------------------------
#SkVent data
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/skvent")

#1) Landmarks
  #land <- readland.tps(file="SkVent_30_10_13_13landmarks+1curve_edited.TPS")
#2) Sliders
  #curves <- as.matrix(read.table(file="SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
#3) Taxonomy
  #taxa <- read.csv("SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
#4) Specimens to remove
  #rem <- read.csv("SkVent_remove_spec.csv", header=T)
#------------------------------------------
#Mandibles data
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/mands")

#1) Landmarks
  #land <- readland.tps(file="Mands_14_03_2014_7landmarks+4curves_edited.TPS")
#2) Sliders
  #curves <- as.matrix(read.table("Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
#3) Taxonomy
  #taxa <- read.csv("Mands_14_03_2014_Images+Specimens.csv", header=T)
#4) Specimens to remove
  #rem <- read.csv("Mands_remove_spec.csv", header=T)

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
#OPTIONS; depending on the data and the analysis
#************************************
#1) Clean up the data
#Clean up the sklat, skvent and mands data
  #doesn't apply to the skdors data because rem is NULL

#find the ID numbers of specimens with missing data
  #matching <- matching.id(rem$SpecID, combine$SpecID)
    #combine <- remove.from.list(combine, matching)
    #combine <- droplevels.from.list(combine)
#*********************************************
 #2) Select which families to work with 
  #2a) Select the tenrec and golden mole specimens only
    tc.gm <- c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))

    
    mydata <- select.from.list(combine, tc.gm)
    mydata <- droplevels.from.list(mydata)
  
  #2b) All of the families (just rename combine as mydata)
    #NB: remove the Notoryctidae because there's only one specimen in the family
      #mydata <- remove.from.list(combine, which(combine$Fam=="Notoryctidae"))
      #mydata <- droplevels.from.list(mydata)
#**************************************************
  #3) Option to remove the Microgale tenrecs

   #mic <- which(mydata$Genus=="Microgale")

   #mydata <- remove.from.list(mydata, mic)
   #mydata <- droplevels.from.list(mydata)
#**************************************************  

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

#PCA graph in the output files section at the end of the script
#######################################
#SELECT PC AXES
#######################################

PC95axes <- selectPCaxes(sps.meanPCA, 0.956, binom)

#select the axes for each family
  gmolePC <- PC95axes[which(sp.fam$Family=="Chrysochloridae"),]
  tenrecPC <- PC95axes[which(sp.fam$Family=="Tenrecidae"),]

#Additional families for the full data set (including Solenodontidae even though there are only two species)
  #hedgePC <- PC95axes[which(sp.fam$Family=="Erinaceidae"),]
  #shrewPC <- PC95axes[which(sp.fam$Family=="Soricidae"),]
  #molePC <- PC95axes[which(sp.fam$Family=="Talpidae"),]
  #solenPC <- PC95axes[which(sp.fam$Family=="Solenodontidae"),]

#######################################
#CALCULATE DISPARITY
#######################################
#Based on PC axes
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

#Based on sum of squared distances (Zeldich 2012)
  #interlandmark distance: compare each species to the overall mean shape of all species
    ild.distance <- dist.to.ref(sps.mean$meanshape)
    #create a matrix version as well with species as rownames (need it for permutation tests later in the script)
    ild.distance.mat <- as.matrix(ild.distance)
    rownames(ild.distance.mat) <- binom
#-----------------------------------------------------    
#Distance to the golden moles meanshape (needed this to send mandible information to Gary Broner) 
  #gmole.sps.mean <- sps.mean$meanshape[,,which(sp.fam$Fam == "Chrysochloridae")]
  #ild.gmole <- dist.to.ref(gmole.sps.mean)
  #ild.gmole.mat <- as.matrix(ild.gmole)
  #rownames(ild.gmole.mat) <- binom[1:12]
#------------------------------------------------------

  #tenrecs
    tenrec.ild <- ild.distance[which(sp.fam$Fam == "Tenrecidae")]
    tenrecMD<- ZelditchMD(tenrec.ild)

  #golden moles
    gmole.ild <- ild.distance[which(sp.fam$Fam == "Chrysochloridae")]
    gmoleMD <- ZelditchMD(gmole.ild)
#------------------------------
#Extra families for the full analysis
  #Hedgehogs
    #variance
    #hedge.v <- PCvariance(hedgePC)
    #hedge.sv <- PCsumvar(hedgePC)
    #hedge.pv <- PCprodvar(hedgePC)

    #range
    #hedge.r <- PCrange(hedgePC)
    #hedge.sr <- PCsumrange(hedgePC)
    #hedge.pr <- PCprodrange(hedgePC)

    #sum of squared Euclidean distances
    #hedge.ild <- ild.distance[which(sp.fam$Fam == "Erinaceidae")]
    #hedgeMD <- ZelditchMD(hedge.ild)             
    
  #Moles
    #variance
    #mole.v <- PCvariance(molePC)
    #mole.sv <- PCsumvar(molePC)
    #mole.pv <- PCprodvar(molePC)

    #range
    #mole.r <- PCrange(molePC)
    #mole.sr <- PCsumrange(molePC)
    #mole.pr <- PCprodrange(molePC)

    #sum of squared Euclidean distances
    #mole.ild <- ild.distance[which(sp.fam$Fam == "Talpidae")]
    #moleMD <- ZelditchMD(mole.ild) 
    
  #Shrews
    #variance
    #shrew.v <- PCvariance(shrewPC)
    #shrew.sv <- PCsumvar(shrewPC)
    #shrew.pv <- PCprodvar(shrewPC)

    #range
    #shrew.r <- PCrange(shrewPC)
    #shrew.sr <- PCsumrange(shrewPC)
    #shrew.pr <- PCprodrange(shrewPC)

    #sum of squared Euclidean distances
    #shrew.ild <- ild.distance[which(sp.fam$Fam == "Soricidae")]
    #shrewMD <- ZelditchMD(shrew.ild) 
    
  #Solenodons
    #variance
    #solen.v <- PCvariance(solenPC)
    #solen.sv <- PCsumvar(solenPC)
    #solen.pv <- PCprodvar(solenPC)

    #range
    #solen.r <- PCrange(solenPC)
    #solen.sr <- PCsumrange(solenPC)
    #solen.pr <- PCprodrange(solenPC)
    
    #sum of squared Euclidean distances
    #solen.ild <- ild.distance[which(sp.fam$Fam == "Solenodontidae")]
    #solenMD <- ZelditchMD(solen.ild) 
#------------------------------------------------    

#Put the disparity calculations into a single table (just tenrecs and golden moles)
  disp <- matrix(NA,nrow=2, ncol=5)
    rownames(disp) <- c("Tenrec","Gmole")
    colnames(disp) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    disp[,1] <- c(tenrec.sv, gmole.sv)
    disp[,2] <- c(tenrec.pv, gmole.pv)
    disp[,3] <- c(tenrec.sr, gmole.sr)
    disp[,4] <- c(tenrec.pr, gmole.pr)
    disp[,5] <- c(tenrecMD, gmoleMD)
    
#Put the disparity calculations into a single table (full table for all families)
  #disp <- matrix(NA,nrow=6, ncol=5)
    #rownames(disp) <- c("Tenrec","Gmole", "Hedge", "Mole", "Shrew", "Solen")
    #colnames(disp) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    #disp[,1] <- c(tenrec.sv, gmole.sv, hedge.sv, mole.sv, shrew.sv, solen.sv)
    #disp[,2] <- c(tenrec.pv, gmole.pv, hedge.pv, mole.pv, shrew.pv, solen.pv)
    #disp[,3] <- c(tenrec.sr, gmole.sr, hedge.sr, mole.sr, shrew.sr, solen.sr)
    #disp[,4] <- c(tenrec.pr, gmole.pr, hedge.pr, mole.pr, shrew.pr, solen.pr)
    #disp[,5] <- c(tenrecMD, gmoleMD, hedgeMD, moleMD, shrewMD, solenMD)
  

#######################################
#COMPARE FAMILIES
#######################################
#Compare morphospace occupation (not comparing disparity metrics directly)
        #Not very meaningful because it just shows that overall families occupy significantly different areas of morphospace
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
  #tests for PC disparity metrics and Zelditch MD
  #Advantage of this method is that it takes differences in sample size into account
  #Pairwise observed differences in disparity among all family groups
  
  obs.diff.sv <- group.pair.diff(group.identity=rownames(disp), group.values=disp[,1])    
  obs.diff.pv <- group.pair.diff(group.identity=rownames(disp), group.values=disp[,2])   
  obs.diff.sr <- group.pair.diff(group.identity=rownames(disp), group.values=disp[,3])   
  obs.diff.pr <- group.pair.diff(group.identity=rownames(disp), group.values=disp[,4])   
  obs.diff.md <- group.pair.diff(group.identity=rownames(disp), group.values=disp[,5])   
  
#Permutation tests for significant differences in disparity between tenrecs and other family groups
  #1) Tenrecs vs golden moles  
      #Permutation tests 
      tc.gm.perm.sv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCsumvar)
      tc.gm.perm.pv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCprodvar)
      tc.gm.perm.sr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCsumrange)
      tc.gm.perm.pr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes, PCprodrange)  
      tc.gm.perm.md <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, ild.distance.mat, ZelditchMD)  
  
      #test for significant differences
          #(NB: Another way to look at the distribution is with a table showing when the permutated values are greater than the observed)
          #tc.gm.obs.dist.sv <- table(tc.gm.perm.sv >= obs.diff.sv[1,3])
      tc.gm.pvalue.sv <- pvalue.dist(tc.gm.perm.sv, obs.diff.sv[1,3])
      tc.gm.pvalue.pv <- pvalue.dist(tc.gm.perm.pv, obs.diff.pv[1,3])
      tc.gm.pvalue.sr <- pvalue.dist(tc.gm.perm.sr, obs.diff.sr[1,3])
      tc.gm.pvalue.pr <- pvalue.dist(tc.gm.perm.pr, obs.diff.pr[1,3])
      tc.gm.pvalue.md <- pvalue.dist(tc.gm.perm.md, obs.diff.md[1,3])

#Summary table of the results for just tenrecs vc. golden moles
      perm.res.summary <- matrix(NA, nrow=5, ncol=7)
        colnames(perm.res.summary) <- c("metric", "obs.tenrec", "obs.gmole", "obs.diff", "perm.min", "perm.max", "pvalue")
        perm.res.summary[,1] <- c("sumvar", "prodvar", "sumrange", "prodrange", "ZelditchMD")
          #tenrec disparity, rounded to 3 significant figures
        perm.res.summary[,2] <- signif(disp[1,],3)
           #gmole disparity, rounded to 3 significant figures
        perm.res.summary[,3] <- signif(disp[2,],3)
         #observed differences, everything rounded to 3 significant digits
        perm.res.summary[,4] <- c(signif(obs.diff.sv[1,3],3), signif(obs.diff.pv[1,3],3), signif(obs.diff.sr[1,3],3), signif(obs.diff.pr[1,3],3), signif(obs.diff.md[1,3],3))
         #minimum permutated values, everything rounded to 3 significant digits
        perm.res.summary[,5] <- c(signif(min(tc.gm.perm.sv),3), signif(min(tc.gm.perm.pv),3), signif(min(tc.gm.perm.sr),3), signif(min(tc.gm.perm.pr),3), signif(min(tc.gm.perm.md),3))
          #maximum permutated values, everything rounded to 3 significant digits
        perm.res.summary[,6] <- c(signif(max(tc.gm.perm.sv),3), signif(max(tc.gm.perm.pv),3), signif(max(tc.gm.perm.sr),3), signif(max(tc.gm.perm.pr),3), signif(max(tc.gm.perm.md),3))        
          #p values from comparing the observed differences to the distribution of permutated values
        perm.res.summary[,7] <- c(tc.gm.pvalue.sv, tc.gm.pvalue.pv, tc.gm.pvalue.sr, tc.gm.pvalue.pr, tc.gm.pvalue.md)

  perm.res.summary <- as.data.frame(perm.res.summary)


#------------------------------------------------------
#Additional comparisons for a full analysis comparing tenrecs to other family groups as well
#-------------------------------------------------------
  #2) Tenrecs vs. hedgehogs
      #Permutation tests 
      tc.hd.perm.sv <- perm.diff.two.groups(1000, "Erinaceidae", "Tenrecidae", sp.fam, PC95axes, PCsumvar)
      tc.hd.perm.pv <- perm.diff.two.groups(1000, "Erinaceidae", "Tenrecidae", sp.fam, PC95axes, PCprodvar)
      tc.hd.perm.sr <- perm.diff.two.groups(1000, "Erinaceidae", "Tenrecidae", sp.fam, PC95axes, PCsumrange)
      tc.hd.perm.pr <- perm.diff.two.groups(1000, "Erinaceidae", "Tenrecidae", sp.fam, PC95axes, PCprodrange)  
      tc.hd.perm.md <- perm.diff.two.groups(1000, "Erinaceidae", "Tenrecidae", sp.fam, ild.distance.mat, ZelditchMD)  
  
      #test for significant differences
      tc.hd.pvalue.sv <- pvalue.dist(tc.hd.perm.sv, obs.diff.sv[2,3])
      tc.hd.pvalue.pv <- pvalue.dist(tc.hd.perm.pv, obs.diff.pv[2,3])
      tc.hd.pvalue.sr <- pvalue.dist(tc.hd.perm.sr, obs.diff.sr[2,3])
      tc.hd.pvalue.pr <- pvalue.dist(tc.hd.perm.pr, obs.diff.pr[2,3])
      tc.hd.pvalue.md <- pvalue.dist(tc.hd.perm.md, obs.diff.md[2,3])

  #3) Tenrecs vs. moles
      #Permutation tests 
      tc.ml.perm.sv <- perm.diff.two.groups(1000, "Talpidae", "Tenrecidae", sp.fam, PC95axes, PCsumvar)
      tc.ml.perm.pv <- perm.diff.two.groups(1000, "Talpidae", "Tenrecidae", sp.fam, PC95axes, PCprodvar)
      tc.ml.perm.sr <- perm.diff.two.groups(1000, "Talpidae", "Tenrecidae", sp.fam, PC95axes, PCsumrange)
      tc.ml.perm.pr <- perm.diff.two.groups(1000, "Talpidae", "Tenrecidae", sp.fam, PC95axes, PCprodrange)  
      tc.ml.perm.md <- perm.diff.two.groups(1000, "Talpidae", "Tenrecidae", sp.fam, ild.distance.mat, ZelditchMD)  
  
      #test for significant differences
      tc.ml.pvalue.sv <- pvalue.dist(tc.ml.perm.sv, obs.diff.sv[3,3])
      tc.ml.pvalue.pv <- pvalue.dist(tc.ml.perm.pv, obs.diff.pv[3,3])
      tc.ml.pvalue.sr <- pvalue.dist(tc.ml.perm.sr, obs.diff.sr[3,3])
      tc.ml.pvalue.pr <- pvalue.dist(tc.ml.perm.pr, obs.diff.pr[3,3])
      tc.ml.pvalue.md <- pvalue.dist(tc.ml.perm.md, obs.diff.md[3,3])

  #4) Tenrecs vs. shrews
      #Permutation tests 
      tc.sh.perm.sv <- perm.diff.two.groups(1000, "Soricidae", "Tenrecidae", sp.fam, PC95axes, PCsumvar)
      tc.sh.perm.pv <- perm.diff.two.groups(1000, "Soricidae", "Tenrecidae", sp.fam, PC95axes, PCprodvar)
      tc.sh.perm.sr <- perm.diff.two.groups(1000, "Soricidae", "Tenrecidae", sp.fam, PC95axes, PCsumrange)
      tc.sh.perm.pr <- perm.diff.two.groups(1000, "Soricidae", "Tenrecidae", sp.fam, PC95axes, PCprodrange)  
      tc.sh.perm.md <- perm.diff.two.groups(1000, "Soricidae", "Tenrecidae", sp.fam, ild.distance.mat, ZelditchMD)  
  
      #test for significant differences
      tc.sh.pvalue.sv <- pvalue.dist(tc.sh.perm.sv, obs.diff.sv[4,3])
      tc.sh.pvalue.pv <- pvalue.dist(tc.sh.perm.pv, obs.diff.pv[4,3])
      tc.sh.pvalue.sr <- pvalue.dist(tc.sh.perm.sr, obs.diff.sr[4,3])
      tc.sh.pvalue.pr <- pvalue.dist(tc.sh.perm.pr, obs.diff.pr[4,3])
      tc.sh.pvalue.md <- pvalue.dist(tc.sh.perm.md, obs.diff.md[4,3])


#Full summary table of the results (make this shorter when it's just a tenrec vs. golden mole analysis)
  perm.res.summary <- matrix(NA, nrow=20, ncol=7)
    colnames(perm.res.summary) <- c("metric", "fam1", "fam2", "obs.diff", "perm.min", "perm.max", "pvalue")
    
  perm.res.summary[,1] <- c(rep("sumvar",4), rep("prodvar", 4), rep("sumrange", 4), rep("prodrange",4), rep("ZelditchMD",4))
  perm.res.summary[,2] <- rep(c("gmole", "hedge", "mole", "shrew"),5)
  perm.res.summary[,3] <- "tenrec"
    #fill in the observed disparity values
      perm.res.summary[1:4,4] <- obs.diff.sv[1:4,3] 
      perm.res.summary[5:8,4] <- obs.diff.pv[1:4,3] 
      perm.res.summary[9:12,4] <- obs.diff.sr[1:4,3] 
      perm.res.summary[13:16,4] <- obs.diff.pr[1:4,3]
      perm.res.summary[17:20,4] <- obs.diff.md[1:4,3]
    #minimum and maximum values from permutations and p values
      #sum of variance
      perm.res.summary[1,5:7] <- c(min(tc.gm.perm.sv), max(tc.gm.perm.sv), tc.gm.pvalue.sv)
      perm.res.summary[2,5:7] <- c(min(tc.hd.perm.sv), max(tc.hd.perm.sv), tc.hd.pvalue.sv)
      perm.res.summary[3,5:7] <- c(min(tc.ml.perm.sv), max(tc.ml.perm.sv), tc.ml.pvalue.sv)
      perm.res.summary[4,5:7] <- c(min(tc.sh.perm.sv), max(tc.sh.perm.sv), tc.sh.pvalue.sv)
      
      #product of variance
      perm.res.summary[5,5:7] <- c(min(tc.gm.perm.pv), max(tc.gm.perm.pv), tc.gm.pvalue.pv)
      perm.res.summary[6,5:7] <- c(min(tc.hd.perm.pv), max(tc.hd.perm.pv), tc.hd.pvalue.pv)
      perm.res.summary[7,5:7] <- c(min(tc.ml.perm.pv), max(tc.ml.perm.pv), tc.ml.pvalue.pv)
      perm.res.summary[8,5:7] <- c(min(tc.sh.perm.pv), max(tc.sh.perm.pv), tc.sh.pvalue.pv)

      #sum of ranges
      perm.res.summary[9,5:7] <- c(min(tc.gm.perm.sr), max(tc.gm.perm.sr), tc.gm.pvalue.sr)
      perm.res.summary[10,5:7] <- c(min(tc.hd.perm.sr), max(tc.hd.perm.sr), tc.hd.pvalue.sr)
      perm.res.summary[11,5:7] <- c(min(tc.ml.perm.sr), max(tc.ml.perm.sr), tc.ml.pvalue.sr)
      perm.res.summary[12,5:7] <- c(min(tc.sh.perm.sr), max(tc.sh.perm.sr), tc.sh.pvalue.sr)
                               
      #product of ranges
      perm.res.summary[13,5:7] <- c(min(tc.gm.perm.pr), max(tc.gm.perm.pr), tc.gm.pvalue.pr)
      perm.res.summary[14,5:7] <- c(min(tc.hd.perm.pr), max(tc.hd.perm.pr), tc.hd.pvalue.pr)
      perm.res.summary[15,5:7] <- c(min(tc.ml.perm.pr), max(tc.ml.perm.pr), tc.ml.pvalue.pr)
      perm.res.summary[16,5:7] <- c(min(tc.sh.perm.pr), max(tc.sh.perm.pr), tc.sh.pvalue.pr)                               
     
      #ZelditchMD
      perm.res.summary[17,5:7] <- c(min(tc.gm.perm.md), max(tc.gm.perm.md), tc.gm.pvalue.md)
      perm.res.summary[18,5:7] <- c(min(tc.hd.perm.md), max(tc.hd.perm.md), tc.hd.pvalue.md)
      perm.res.summary[19,5:7] <- c(min(tc.ml.perm.md), max(tc.ml.perm.md), tc.ml.pvalue.md)
      perm.res.summary[20,5:7] <- c(min(tc.sh.perm.md), max(tc.sh.perm.md), tc.sh.pvalue.md)
                               
perm.res.summary<-as.data.frame(perm.res.summary)                                    
  

#Extra comparison of moles vs. golden moles (mainly relevant when I move on to trying to figure out why gmoles have more variable mandibles than tenrecs)
      #Permutation tests 
      gm.ml.perm.sv <- perm.diff.two.groups(1000, "Talpidae", "Chrysochloridae", sp.fam, PC95axes, PCsumvar)
      gm.ml.perm.pv <- perm.diff.two.groups(1000, "Talpidae", "Chrysochloridae", sp.fam, PC95axes, PCprodvar)
      gm.ml.perm.sr <- perm.diff.two.groups(1000, "Talpidae", "Chrysochloridae", sp.fam, PC95axes, PCsumrange)
      gm.ml.perm.pr <- perm.diff.two.groups(1000, "Talpidae", "Chrysochloridae", sp.fam, PC95axes, PCprodrange)  
      gm.ml.perm.md <- perm.diff.two.groups(1000, "Talpidae", "Chrysochloridae", sp.fam, ild.distance.mat, ZelditchMD)  
  
      #test for significant differences
      gm.ml.pvalue.sv <- pvalue.dist(gm.ml.perm.sv, obs.diff.sv[7,3])
      gm.ml.pvalue.pv <- pvalue.dist(gm.ml.perm.pv, obs.diff.pv[7,3])
      gm.ml.pvalue.sr <- pvalue.dist(gm.ml.perm.sr, obs.diff.sr[7,3])
      gm.ml.pvalue.pr <- pvalue.dist(gm.ml.perm.pr, obs.diff.pr[7,3])
      gm.ml.pvalue.md <- pvalue.dist(gm.ml.perm.md, obs.diff.md[7,3])
  #Only the ZelditchMD values are significantly different


#######################################
#OUTPUT FILES
#######################################

#Save the outputs to different working directory
#SkDors
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skdors")
#SkLat
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/sklat")
#SkVent
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skvent")
#Mands
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")

#**************************************************
#1) Shape data, taxonomy , disparity measures, disparity comparisons 
#****************************************************
#SkDors: 
#A) All tenrecs and golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="SkDors_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkDors_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="SkDors_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="SkDors_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
      
      
      
#B) Non-microgale tenrecs and golden moles

  #1) Average shape coordinates
    #dput(sps.mean, file="SkDors_nonmic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkDors_nonmic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="SkDors_nonmic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="SkDors_nonmic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)


#--------------------------------------------------------------------
#SkLat

#A) All tenrecs and golden moles
  #1) Average shape coordinates
     #dput(sps.mean, file="SkLat_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
     #write.table(file="SkLat_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
     #write.table(file="SkLat_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
      #write.table(file="SkLat_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
           
#B) Non-microgale tenrecs and golden moles
  #1) Average shape coordinates
     #dput(sps.mean, file="SkLat_nonmic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
     #write.table(file="SkLat_nonmic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
     #write.table(file="SkLat_nonmic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
     #write.table(file="SkLat_nonmic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)


#----------------------------------------------------------
#SkVent
#A) All tenrecs and golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="SkVent_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkVent_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="SkVent_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="SkVent_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
    
    
#B) Non-microgale tenrecs and golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="SkVent_nonmic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkVent_nonmic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="SkVent_nonmic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="SkVent_nonmic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#----------------------------------------------------------
#Mands
#A) All tenrecs and golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="Mands_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="Mands_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="Mands_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="Mands_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  

#B) Non-microgale tenrecs and golden moles
  #1) Average shape coordinates
    #dput(sps.mean, file="Mands_nonmic_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="Mands_nonmic_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of npMANOVA results
    #write.table(file="Mands_nonmic_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of permutation test for significant differences in group disparity
    #write.table(file="Mands_nonmic_tenrec+gmole_disp.signif.txt",perm.res.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#**************************************************
#2) Save the PCA plot
#******************************************
  #SkDors
    #pdf(file="skdors_tenrec+gmole_PCA.pdf")
    #pdf(file="skdors_nonmic_tenrec+gmole_PCA.pdf")
  
  #SkLat
    #pdf(file="sklat_tenrec+gmole_PCA.pdf") 
    #pdf(file="sklat_nonmic_tenrec+gmole_PCA.pdf") 
  
  #SkVent
    #pdf(file="skvent_tenrec+gmole_PCA.pdf")
    #pdf(file="skvent_nonmic_tenrec+gmole_PCA.pdf")
        
  #Mands
   #pdf(file="mands_tenrec+gmole_PCA.pdf")
   #pdf(file="mands_nonmic_tenrec+gmole_PCA.pdf")

#PCA plot, default colour palette so Chrysochloridae are black and Tenrecidae are red
 
  #plot(xaxis,yaxis, xlab="Species' average PC1", ylab="Species' average PC2",las=1,
       #col=sp.fam$Family,pch=16, bty="l",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
    #draw the min,max and 0 values on the x axis
      #axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
    #same for the y axis
      #axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
    #add dotted lines along 0,0
      #abline(0,0,h=0,v=0,lty=2,lwd=1)
      
#dev.off()

#identify points on the graph
#  identify(xaxis,yaxis,labels=(sp.fam$Binom))

#------------------------------------------------------  
#Cleaned up plot for presentations (no axes or values)
  #Plot just the golden moles on their own (tenrecs are white)
#  gmole.alone <- c("black", "white")
#  palette(gmole.alone)

#  plot(xaxis,yaxis, xlab="", ylab="",las=1,
#       col=sp.fam$Family,pch=16, bty="n",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
#      abline(0,0,h=0,v=0,lty=1,lwd=2)

  #Plot golden moles and tenrecs
#  palette("default")

#  plot(xaxis,yaxis, xlab="", ylab="",las=1,
#       col=sp.fam$Family,pch=16, bty="n",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
#      abline(0,0,h=0,v=0,lty=1,lwd=2)





