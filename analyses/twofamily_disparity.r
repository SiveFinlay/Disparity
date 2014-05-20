#09/05/2014

#Re-coding and re-analysis of my disparity project

#general script which can be used with any of the data sets (skdors, skvent, sklat, mands)
  #change the input files depending on which data I'm using
  
#compare disparity in tenrecs and golden moles only
  #additional option of selecting just microgale tenrecs

#steps:
  #1) Read in a clean up raw landmark data to select tenrecs and golden moles only
  #2) Procrustes superimposition of tenrecs and golden moles
  #3) Find the average Procrustes shape coordinates for each species
  #4) PCA of the shape coordinates
  #5) Select PC axes that account for 95% of the variation
  #6) Calculate disparity measures
  #7) Compare disparity in families; npMANOVA
  #8) Output files: shape data, taxonomy, disparity, disparity comparisons
  #8) Sensitivity analysis (rarefaction)

#output from this script
  #data
    #1) The average shape coordinates for each species of the GPA-aligned specimens (sps.mean)
    #2) The taxonomic information (Family and Binomial) for these shape coordinates (sp.fam)
        # (the binomial names are in the same order in each object)
    #3) Table of disparity measures of each family                                  (disp)
    #4) Table of npMANOVA reults; based on distance matrix and PC axes              (manova.res)
    
    
  #figures
    #1) PCA plots
    #2) Rarefaction profiles
    
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
  #setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/skdors")

#1) Landmarks
#landmarks + curves file with the control lines removed
  #land <- readland.tps(file="Skdors_16_12_13_10landmarks+4curves_edited.TPS")

#2) Sliders
#edited sliders file (top 2 rows removed and the words before slide after put in instead
  #curves <- as.matrix(read.table("Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
#file that has the correct taxonomy for each of the images
  #taxa <- read.csv ("Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)

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
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/mands")

#1) Landmarks
  land <- readland.tps(file="Mands_14_03_2014_7landmarks+4curves_edited.TPS")
#2) Sliders
  curves <- as.matrix(read.table("Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
#3) Taxonomy
  taxa <- read.csv("Mands_14_03_2014_Images+Specimens.csv", header=T)
#4) Specimens to remove
  rem <- read.csv("Mands_remove_spec.csv", header=T)

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
#Option depending on the data
#Clean up the sklat, skvent and mands data
  #doesn't apply to the skdors data because rem is NULL
#************************************************** 
#find the ID numbers of specimens with missing data
  matching <- matching.id(rem$SpecID, combine$SpecID)
    combine <- remove.from.list(combine, matching)
    combine <- droplevels.from.list(combine)
#*********************************************

#Select the tenrec and golden mole specimens only
  tc.gm <- c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))
  
  mydata <- select.from.list(combine, tc.gm)
  mydata <- droplevels.from.list(mydata)

#**************************************************
#Option depending on the analysis
#************************************************** 
#Option to remove all of the Microgale specimens
  # mic <- which(mydata$Genus=="Microgale")

  # mydata <- remove.from.list(mydata, mic)
  # mydata <- droplevels.from.list(mydata)
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
  xaxis <- sps.meanPCA$pc.scores[,1]
  yaxis <- sps.meanPCA$pc.scores[,2]

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

  #tenrecs
    tenrec.ild <- ild.distance[which(sp.fam$Fam == "Tenrecidae")]
    tenrecMD<- ZelditchMD(tenrec.ild)

  #golden moles
    gmole.ild <- ild.distance[which(sp.fam$Fam == "Chrysochloridae")]
    gmoleMD <- ZelditchMD(gmole.ild)

#Put the disparity calculations into a single table
  disp <- matrix(NA,nrow=2, ncol=5)
    rownames(disp) <- c("Tenrec","Gmole")
    colnames(disp) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    disp[,1] <- c(tenrec.sv,gmole.sv)
    disp[,2] <- c(tenrec.pv,gmole.pv)
    disp[,3] <- c(tenrec.sr,gmole.sr)
    disp[,4] <- c(tenrec.pr,gmole.pr)
    disp[,5] <- c(tenrecMD,gmoleMD)
  

#######################################
#COMPARE FAMILIES
#######################################
#Compare morphospace occupation (not comparing disparity metrics directly)

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
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")

#**************************************************
#1) Shape data, taxonomy , disparity measures, disparity comparisons 
#SkDors
  #1) Average shape coordinates
    #dput(sps.mean, file="SkDors_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkDors_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of disparity measures for each family
    #write.table(file="SkDors_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of npMANOVA results
    #write.table(file="SkDors_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#--------------------------------------------------------------------
#SkLat
  #1) Average shape coordinates
     #dput(sps.mean, file="SkLat_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
     #write.table(file="SkLat_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of disparity measures for each family
     #write.table(file="SkLat_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of npMANOVA results
     #write.table(file="SkLat_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#----------------------------------------------------------
#SkVent
  #1) Average shape coordinates
    #dput(sps.mean, file="SkVent_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    #write.table(file="SkVent_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of disparity measures for each family
    #write.table(file="SkVent_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of npMANOVA results
    #write.table(file="SkVent_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
#----------------------------------------------------------
#Mands
  #1) Average shape coordinates
#    dput(sps.mean, file="Mands_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
#    write.table(file="Mands_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #3) Table of disparity measures for each family
#    write.table(file="Mands_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  #4) Table of npMANOVA results
#    write.table(file="Mands_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)

#**************************************************
#2) Save the PCA plot
  #SkDors
    #pdf(file="skdors_tenrec+gmole_PCA.pdf")
  #SkLat
    #pdf(file="sklat_tenrec+gmole_PCA.pdf") 
  #SkVent
    #pdf(file="skvent_tenrec+gmole_PCA.pdf")
  #Mands
    pdf(file="mands_tenrec+gmole_PCA.pdf")

#PCA plot, default colour palette so Chrysochloridae are black and Tenrecidae are red
 
  plot(xaxis,yaxis, xlab="Species' average PC1", ylab="Species' average PC2",las=1,
       col=sp.fam$Family,pch=16, bty="l",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
    #same for the y axis
      axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
    #add dotted lines along 0,0
      abline(0,0,h=0,v=0,lty=2,lwd=1)

dev.off()
#identify points on the graph
  #identify(xaxis,yaxis,labels=(sp.fam$Binom))

#########################################################################
#SENSITIVITY ANALYSIS for the PC-based disparity metrics
#########################################################################

#Re-sample the PC axis data: from 2 to (full species-1), resample without replacement 
#Tenrecs 
  tenrec.res <- resample.data(tenrecPC, samp.min=2, samp.max=(nrow(tenrecPC)-1), no.replicates =100, no.col=ncol(tenrecPC))

#Golden moles
  gmole.res <- resample.data(gmolePC, samp.min=2, samp.max=(nrow(gmolePC)-1), no.replicates =100, no.col=ncol(gmolePC))

#-----------------------------------------------------------------
#Calculate disparity metrics for each re-sampled data
  #Sum of variance
    tenrec.res.sv <- calc.each.array(tenrec.res, PCsumvar)
    gmole.res.sv <- calc.each.array(gmole.res, PCsumvar)
  
  #Product of variance
    tenrec.res.pv <- calc.each.array(tenrec.res, PCprodvar)
    gmole.res.pv <- calc.each.array(gmole.res, PCprodvar)
  
  #Sum of ranges
    tenrec.res.sr <- calc.each.array(tenrec.res, PCsumrange)
    gmole.res.sr <- calc.each.array(gmole.res, PCsumrange)
  
  #Product of ranges
    tenrec.res.pr <- calc.each.array(tenrec.res, PCprodrange)
    gmole.res.pr <- calc.each.array(gmole.res, PCprodrange)

#-----------------------------------------------------------------
#Mean values for each sample (mean sumvar for 2 species, 3 species etc.)
  #Sum of variance
    tenrec.res.sv.mean <- lapply(tenrec.res.sv, mean)
    gmole.res.sv.mean <- lapply(gmole.res.sv, mean)

  #Product of variance
    tenrec.res.pv.mean <- lapply(tenrec.res.pv, mean)
    gmole.res.pv.mean <- lapply(gmole.res.pv, mean)
    
  #Sum of ranges
    tenrec.res.sr.mean <- lapply(tenrec.res.sr, mean)
    gmole.res.sr.mean <- lapply(gmole.res.sr, mean)
    
  #Product of ranges
    tenrec.res.pr.mean <- lapply(tenrec.res.pr, mean)
    gmole.res.pr.mean <- lapply(gmole.res.pr, mean)

#----------------------------------------------------
#Bootstrap confidence intervals
  #Sum of variance
    #Tenrecs
      tenrec.sv.boot <- lapply(tenrec.res.sv, boot.mean.1000)
      tenrec.sv.min.conf <- unlist(lapply(tenrec.sv.boot, boot.95.min.confidence))
      tenrec.sv.max.conf <- unlist(lapply(tenrec.sv.boot, boot.95.max.confidence))
    #Golden moles
      gmole.sv.boot <- lapply(gmole.res.sv, boot.mean.1000)
      gmole.sv.min.conf <- unlist(lapply(gmole.sv.boot, boot.95.min.confidence))
      gmole.sv.max.conf <- unlist(lapply(gmole.sv.boot, boot.95.max.confidence))      

  #Product of variance
    #Tenrecs
      tenrec.pv.boot <- lapply(tenrec.res.pv, boot.mean.1000)
      tenrec.pv.min.conf <- unlist(lapply(tenrec.pv.boot, boot.95.min.confidence))
      tenrec.pv.max.conf <- unlist(lapply(tenrec.pv.boot, boot.95.max.confidence))
    #Golden moles
      gmole.pv.boot <- lapply(gmole.res.pv, boot.mean.1000)
      gmole.pv.min.conf <- unlist(lapply(gmole.pv.boot, boot.95.min.confidence))
      gmole.pv.max.conf <- unlist(lapply(gmole.pv.boot, boot.95.max.confidence))      

  #Sum of ranges
    #Tenrecs
      tenrec.sr.boot <- lapply(tenrec.res.sr, boot.mean.1000)
      tenrec.sr.min.conf <- unlist(lapply(tenrec.sr.boot, boot.95.min.confidence))
      tenrec.sr.max.conf <- unlist(lapply(tenrec.sr.boot, boot.95.max.confidence))
    #Golden moles
      gmole.sr.boot <- lapply(gmole.res.sr, boot.mean.1000)
      gmole.sr.min.conf <- unlist(lapply(gmole.sr.boot, boot.95.min.confidence))
      gmole.sr.max.conf <- unlist(lapply(gmole.sr.boot, boot.95.max.confidence))      

  #Product of ranges
    #Tenrecs
      tenrec.pr.boot <- lapply(tenrec.res.pr, boot.mean.1000)
      tenrec.pr.min.conf <- unlist(lapply(tenrec.pr.boot, boot.95.min.confidence))
      tenrec.pr.max.conf <- unlist(lapply(tenrec.pr.boot, boot.95.max.confidence))
    #Golden moles
      gmole.pr.boot <- lapply(gmole.res.pr, boot.mean.1000)
      gmole.pr.min.conf <- unlist(lapply(gmole.pr.boot, boot.95.min.confidence))
      gmole.pr.max.conf <- unlist(lapply(gmole.pr.boot, boot.95.max.confidence))      

#-------------------------------------------------------------
#Plot the rarefaction curves
  #sample sizes
  tenrec.samp <- c(2:(nrow(tenrecPC)-1)) 
  gmole.samp <-  c(2:(nrow(gmolePC)-1))  

#Range of confidence intervals to use as the ylim values
  conf.range.sv <- c(min(gmole.sv.min.conf), max(gmole.sv.max.conf), min(tenrec.sv.min.conf), max(tenrec.sv.max.conf))
  conf.range.pv <- c(min(gmole.pv.min.conf), max(gmole.pv.max.conf), min(tenrec.pv.min.conf), max(tenrec.pv.max.conf))
  conf.range.sr <- c(min(gmole.sr.min.conf), max(gmole.sr.max.conf), min(tenrec.sr.min.conf), max(tenrec.sr.max.conf))
  conf.range.pr <- c(min(gmole.pr.min.conf), max(gmole.pr.max.conf), min(tenrec.pr.min.conf), max(tenrec.pr.max.conf))
  
  
#SkDors
   #pdf(file="skdors_trc+gmole_PCrarefaction.pdf")

#SkLat
   #pdf(file="sklat_trc+gmole_PCrarefaction.pdf")

#SkVent
   #pdf(file="skvent_trc+gmole_PCrarefaction.pdf")

#Mands
   pdf(file="mands_trc+gmole_PCrarefaction.pdf")


par(mfrow=c(2,2))
#Sum of variance
  plot(tenrec.samp,tenrec.res.sv.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.65, cex.lab=0.8,
       ylim=c(min(conf.range.sv), max(conf.range.sv)),  #range from the lowest minimum confidence interval to the highest maximum confidence interval      
       xlab="SampleSize", ylab="Sum of Variance")
        
         lines(tenrec.samp,tenrec.sv.min.conf, type="o", col="pink")   #confidence intervals for tenrecs
         lines(tenrec.samp,tenrec.sv.max.conf, type="o", col="pink")
   
       lines(gmole.samp,gmole.res.sv.mean,type="o", col="black")  #mean line for golden moles
         lines(gmole.samp,gmole.sv.min.conf, type="o", col="grey") #confidence intervals for golden moles
         lines(gmole.samp,gmole.sv.max.conf, type="o", col="grey")
            
#Product of variance
  plot(tenrec.samp,tenrec.res.pv.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.65, cex.lab=0.8,
       ylim=c(min(conf.range.pv), max(conf.range.pv)),        
       xlab="SampleSize", ylab="Product of Variance")
        
         lines(tenrec.samp,tenrec.pv.min.conf, type="o", col="pink")   
         lines(tenrec.samp,tenrec.pv.max.conf, type="o", col="pink")
   
        lines(gmole.samp,gmole.res.pv.mean,type="o", col="black")  
          lines(gmole.samp,gmole.pv.min.conf, type="o", col="grey") 
          lines(gmole.samp,gmole.pv.max.conf, type="o", col="grey")
            
#Sum of ranges
  plot(tenrec.samp,tenrec.res.sr.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.65, cex.lab=0.8,
       ylim=c(min(conf.range.sr), max(conf.range.sr)),       
       xlab="SampleSize", ylab="Sum of Ranges")
        
         lines(tenrec.samp,tenrec.sr.min.conf, type="o", col="pink")   
         lines(tenrec.samp,tenrec.sr.max.conf, type="o", col="pink")
   
       lines(gmole.samp,gmole.res.sr.mean,type="o", col="black")  
         lines(gmole.samp,gmole.sr.min.conf, type="o", col="grey") 
         lines(gmole.samp,gmole.sr.max.conf, type="o", col="grey")
            
#Product of ranges
  plot(tenrec.samp,tenrec.res.pr.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.65, cex.lab=0.8,
       ylim=c(min(conf.range.pr), max(conf.range.pr)),        
       xlab="SampleSize", ylab="Product of Ranges")
        
         lines(tenrec.samp,tenrec.pr.min.conf, type="o", col="pink")   
         lines(tenrec.samp,tenrec.pr.max.conf, type="o", col="pink")
   
        lines(gmole.samp,gmole.res.pr.mean,type="o", col="black")  
          lines(gmole.samp,gmole.pr.min.conf, type="o", col="grey") 
          lines(gmole.samp,gmole.pr.max.conf, type="o", col="grey")
            
dev.off()    


########################################################
#Come back to this

#Resample the tenrec data with replacment
  #resample with replacement
  tenrec.res.rep <- resample.data(tenrecPC, samp.min=2, samp.max=(nrow(tenrecPC)-1), no.replicates=100, no.col=ncol(tenrecPC))
  #Sum of variance
    tenrec.res.rep.sv <- calc.each.array(tenrec.res.rep, PCsumvar)
    #mean values Sum of variance
    tenrec.res.rep.sv.mean <- lapply(tenrec.res.rep.sv, mean)
    #bootstrap confidence intervals
    #Tenrecs
      tenrec.rep.sv.boot <- lapply(tenrec.res.rep.sv, boot.mean.1000)
      tenrec.rep.sv.min.conf <- unlist(lapply(tenrec.rep.sv.boot, boot.95.min.confidence))
      tenrec.rep.sv.max.conf <- unlist(lapply(tenrec.rep.sv.boot, boot.95.max.confidence))
      
  plot(tenrec.samp,tenrec.res.rep.sv.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.65, cex.lab=0.8,
       ylim=c(min(tenrec.rep.sv.min.conf), max(tenrec.rep.sv.max.conf)),  #range from the lowest minimum confidence interval to the highest maximum confidence interval      
       xlab="SampleSize", ylab="Sum of Variance")
        
         lines(tenrec.samp,tenrec.rep.sv.min.conf, type="o", col="pink")   #confidence intervals for tenrecs
         lines(tenrec.samp,tenrec.rep.sv.max.conf, type="o", col="pink")

#Alternative way of getting confidence intervals (Foote 1992)
#Trial way to get 90% confidence intervals (Foote 1992)
   
#   sorted.tc.sv <- NULL
#    for (i in 1:length(tenrec.res.sv)){
#      sorted.tc.sv[[i]] <- sort(tenrec.res.sv[[i]])
#    }
    
#    min.lim.tc.sv <- NULL
#      for (i in 1:length(sorted.tc.sv)){
#        min.lim.tc.sv[[i]] <- sorted.tc.sv[[i]][6]
#      }
      
#    max.lim.tc.sv <- NULL
#      for (i in 1:length(sorted.tc.sv)){
#        max.lim.tc.sv[[i]] <- sorted.tc.sv[[i]][95]
#      }
#---------------------------    
#    sorted.gm.sv <- NULL
#     for (i in 1:length(gmole.res.sv)){
#      sorted.gm.sv[[i]] <- sort(gmole.res.sv[[i]])
#    }  

#  min.lim.gm.sv <- NULL
#      for (i in 1:length(sorted.gm.sv)){
#        min.lim.gm.sv[[i]] <- sorted.gm.sv[[i]][6]
#      }
      
#    max.lim.gm.sv <- NULL
#      for (i in 1:length(sorted.gm.sv)){
#        max.lim.gm.sv[[i]] <- sorted.gm.sv[[i]][95]
#      }

#tc.samp <- 2:30
#gm.samp <- 2:11 
#dev.new()
#plot(tc.samp,tenrec.res.sv.mean,type="o", bty ="l", las=1,col="red", cex.axis=0.75,
#      ylim=c((min(min.lim.tc.sv)), (max(max.lim.tc.sv))),  #range from the lowest minimum confidence interval to the highest maximum confidence interval      
#      xlab="SampleSize", ylab="Sum of Variance")
#        lines(tc.samp,max.lim.tc.sv, type="o", col="pink")   #confidence intervals for tenrecs
#        lines(tc.samp,min.lim.tc.sv, type="o", col="pink")
   
#    lines(gm.samp,gmole.res.sv.mean,type="o", col="black")  #mean line for golden moles
#            lines(gm.samp,max.lim.gm.sv, type="o", col="grey") #confidence intervals for golden moles
#            lines(gm.samp,min.lim.gm.sv, type="o", col="grey")
            

