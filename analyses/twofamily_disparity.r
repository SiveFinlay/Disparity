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
#----------------------------------------------------------------------------------
#Read in data; directory will change for each data set
#-------------------------------------------------------------
#SkDors data
#    setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/skdors")

#1) Landmarks
#landmarks + curves file with the control lines removed
#   land<-readland.tps(file="Skdors_16_12_13_10landmarks+4curves_edited.TPS")

#2) Sliders
#edited sliders file (top 2 rows removed and the words before slide after put in instead
#    curves<-as.matrix(read.table("Skdors_16_12_13_10landmarks+4curves_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
#file that has the correct taxonomy for each of the images
#    taxa<-read.csv ("Skdors_16_12_13_10landmarks_images+specimens.csv" , header=T)

#4) Specimens to remove
#    Null
#--------------------------------------------------------
#SkLat data
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/sklat")

#1) Landmarks
#    land<-readland.tps(file="SkLat_08_11_13_9landmarks_2curves_edited.TPS")
#2) Sliders
#   curves<-as.matrix(read.table(file="SkLat_08_11_13_9landmarks_2curves_sliders_edited.NTS", header=TRUE))
#3) Taxonomy
#   taxa<-read.csv("SkLat_08_11_13_Specimens+images.csv", header=TRUE)
#4) Specimens to remove
#   rem<-read.csv("SkLat_remove_spec.csv", header=T)

#-----------------------------------------------------
#SkVent data
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/skvent")

#1) Landmarks
#   land<-readland.tps(file="SkVent_30_10_13_13landmarks+1curve_edited.TPS")
#2) Sliders
#   curves<-as.matrix(read.table(file="SkVent_1skull_13landmarks+1curve_sliders_edited.tps", header=TRUE))     #this is a tps file and the others are nts but it doesn't make a difference
#3) Taxonomy
#   taxa<-read.csv("SkVent_30_10_13_imagelist+specimens.csv" , header=TRUE)
#4) Specimens to remove
#   rem<-read.csv("SkVent_remove_spec.csv", header=T)
#------------------------------------------
#Mandibles data
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/mands")

#1) Landmarks
   land<-readland.tps(file="Mands_14_03_2014_7landmarks+4curves_edited.TPS")
#2) Sliders
   curves<-as.matrix(read.table("Mands_14_03_2014_7landmarks+4curves_sliders_edited.txt", header=TRUE))
#3) Taxonomy
   taxa<-read.csv("Mands_14_03_2014_Images+Specimens.csv", header=T)
#4) Specimens to remove
   rem<-read.csv("Mands_remove_spec.csv", header=T)

###############################################################################
#CLEAN UP THE DATA
#################################################
#Combine the taxonomic information with the array of landmark data

combine<-list(land=land, curves=curves, ID=taxa$ID,SpecID=taxa$SpecID, Order=taxa$Order_05, Fam=taxa$Family_05,
          Genus=taxa$Genus_05, Species=taxa$Species_05, Binom=taxa$Binomial_05)

# Remove the _sp. specimens
sp<-which(combine$Species=="sp.")

combine<-remove.from.list(combine, sp)
combine<-droplevels.from.list(combine)

#******************************************************
#Option depending on the data
#Clean up the sklat, skvent and mands data
#*********************************************
#find the ID numbers of specimens with missing data
  #doesn't apply to the skdors data because rem is NULL
  
  matching<-matching.id(rem$SpecID, combine$SpecID)

   combine<-remove.from.list(combine, matching)
   combine<-droplevels.from.list(combine)

#-------------------------------------------------------

#Select the tenrec and golden mole specimens only

tc.gm<-c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))

mydata<-select.from.list(combine, tc.gm)
mydata<-droplevels.from.list(mydata)

#*******************************************************
#Option depending on the analysis
#**************************************************
#Option to remove all of the Microgale specimens
  # mic<-which(mydata$Genus=="Microgale")

  # mydata<-remove.from.list(mydata, mic)
  # mydata<-droplevels.from.list(mydata)
  
#######################################
#PROCRUSTES SUPERIMPOSTION
#######################################

#General Procrustes Alignment of all of the scaled coordinates
mydataGPA<-gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE,)
  #ProcD=TRUE means that the coordinates are aligned by procrustes distance rather than bending energy
      # which means that RWA is equivalent to PCA (Zelditch 2012, page 150)

#List the coordinates with the taxonomic information
Proc.co<-list(coords=mydataGPA$coords,csize=mydataGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
              Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

#######################################
#SPECIES AVERAGING
#######################################

#group the arrays of coordinates according to species
group.sps.coords<-species.coordinates(Proc.co$coords, Proc.co$Binom)

#average coordinate values for each species
sps.mean<-mean.coords(group.sps.coords)

#list of species
binom<-sps.mean$Binom
#######################################
#PRINCIPAL COMPONENTS ANALYSIS
#######################################

sps.meanPCA<-plotTangentSpace(sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

#--------------------------------------------------
#Re-create the PCA plot
xaxis<-sps.meanPCA$pc.scores[,1]
yaxis<-sps.meanPCA$pc.scores[,2]

#colour points by family
  #data frame of the unique family and binomical combinations
  sp.fam<-as.data.frame(unique(cbind(as.matrix(Proc.co$Fam), as.matrix(Proc.co$Binom))))
  colnames(sp.fam)<-c("Family","Binomial")

#PCA plot, default colour palette so Chrysochloridae are black and Tenrecidae are red
dev.new()
plot(xaxis,yaxis, xlab="Species' average PC1", ylab="Species' average PC2",las=1,
  col=sp.fam$Family,pch=16, bty="l",cex.lab=1.5,cex=1.2, xaxt="n",yaxt="n")
    #draw the min,max and 0 values on the x axis
    axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
    #same for the y axis
    axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
  #add dotted lines along 0,0
  abline(0,0,h=0,v=0,lty=2,lwd=1.5)

#identify points on the graph
#identify(xaxis,yaxis,labels=(sp.fam$Binom))

#######################################
#SELECT PC AXES
#######################################

PC95axes<-selectPCaxes(sps.meanPCA, 0.956, binom)

#select the axes for each family
gmolePC<-PC95axes[which(sp.fam$Family=="Chrysochloridae"),]
tenrecPC<-PC95axes[which(sp.fam$Family=="Tenrecidae"),]

#######################################
#CALCULATE DISPARITY
#######################################
#Based on PC axes
  #Tenrecs
  #variance
  tenrec.v<-PCvariance(tenrecPC)
  tenrec.sv<-PCsumvar(tenrecPC)
  tenrec.pv<-PCprodvar(tenrecPC)

  #range
  tenrec.r<-PCrange(tenrecPC)
  tenrec.sr<-PCsumrange(tenrecPC)
  tenrec.pr<-PCprodrange(tenrecPC)


#Disparity measures for the golden moles
  #variance
  gmole.v<-PCvariance(gmolePC)
  gmole.sv<-PCsumvar(gmolePC)
  gmole.pv<-PCprodvar(gmolePC)

  #range
  gmole.r<-PCrange(gmolePC)
  gmole.sr<-PCsumrange(gmolePC)
  gmole.pr<-PCprodrange(gmolePC)

#Based on sum of squared distances (Zeldich 2012)

#interlandmark distance: compare each species to the overall mean shape of all species
ild.distance<-dist.to.ref(sps.mean$meanshape, sp.fam$Fam, sp.fam$Binom)

  #tenrecs
    tenrec.ild<-subset.matrix(ild.distance, ild.distance$Fam, "Tenrecidae")
    tenrecMD<-ZelditchMD(tenrec.ild$Ild)

  #golden moles
    gmole.ild<-subset.matrix(ild.distance, ild.distance$Fam, "Chrysochloridae")
    gmoleMD<-ZelditchMD(gmole.ild$Ild)

#Put the disparity calculations into a single table
disp<-matrix(NA,nrow=2, ncol=5)
  rownames(disp)<-c("Tenrec","Gmole")
  colnames(disp)<-c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
  disp[,1]<-c(tenrec.sv,gmole.sv)
  disp[,2]<-c(tenrec.pv,gmole.pv)
  disp[,3]<-c(tenrec.sr,gmole.sr)
  disp[,4]<-c(tenrec.pr,gmole.pr)
  disp[,5]<-c(tenrecMD,gmoleMD)
  

#######################################
#COMPARE FAMILIES
#######################################
#Compare morphospace occupation (not comparing disparity metrics directly)

#1) Distance matrix

#Euclidean distance matrix of all of the species
  Euc.dist<-as.matrix(dist(PC95axes,method="euclidean", diag=FALSE,upper=FALSE))


#npMANOVA of the distance matrix separated by family
  dist.man<-adonis(Euc.dist~sp.fam$Family, data=sp.fam, permutations=9999,method="euclidean")

#extract the f, r2 and p values
  dist.man.frp<-anova.frp(dist.man)

#NPMANOVA of the PC axes  (e.g. Stayton 2005 and Ruta 2013)
  PC.man<-adonis(PC95axes~sp.fam$Family, data=sp.fam, permutations=999, method="euclidean")

#extract the f, r2 and p values
  PC.man.frp<-anova.frp(PC.man)

#Compare the distance and PC npMANOVA results
  manova.res<-rbind(dist.man.frp, PC.man.frp)
  rownames(manova.res)<-c("dist.man", "PC.man")
#######################################
#SENSITIVITY ANALYSIS
#######################################


#######################################
#OUTPUT FILES
#######################################

#Save the outputs to different working directory

#SkDors
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skdors")
  #1) Average shape coordinates
#    dput(sps.mean, file="SkDors_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
#    write.table(file="SkDors_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #3) Table of disparity measures for each family
#    write.table(file="SkDors_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #4) Table of npMANOVA results
#    write.table(file="SkDors_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#--------------------------------------------------------------------
#SkLat
 # setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/sklat")
  #1) Average shape coordinates
 #   dput(sps.mean, file="SkLat_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
 #   write.table(file="SkLat_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #3) Table of disparity measures for each family
 #   write.table(file="SkLat_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #4) Table of npMANOVA results
 #   write.table(file="SkLat_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#----------------------------------------------------------
#SkVent
#  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/skvent")
  #1) Average shape coordinates
#    dput(sps.mean, file="SkVent_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
#    write.table(file="SkVent_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #3) Table of disparity measures for each family
#    write.table(file="SkVent_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #4) Table of npMANOVA results
#    write.table(file="SkVent_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=T)
#----------------------------------------------------------
#Mands
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")
  #1) Average shape coordinates
    dput(sps.mean, file="Mands_tenrec+gmole_sps.mean.txt")
  #2) Family and species taxonomy
    write.table(file="Mands_tenrec+gmole_sps.mean_taxonomy.txt",sp.fam,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #3) Table of disparity measures for each family
    write.table(file="Mands_tenrec+gmole_disp.txt",disp,col.names=T, row.names=T,sep="\t",quote=F,append=T)
  #4) Table of npMANOVA results
    write.table(file="Mands_tenrec+gmole_manova.res.txt",manova.res,col.names=T, row.names=T,sep="\t",quote=F,append=T)

