#17/07/2014
#Different analyses of disparity in just the mandibles
  #Full analysis showed greater diversity in golden moles than tenrecs
  #Most of the landmarks concentrate around the back(muscle attachement) end of the mandible
    #So the greater diversity might be due to more variable muscle attachment areas in golden moles (not sure why)
    #This script is to test whether you get a different pattern when just looking at the landmarks
      #(NB: Fewer landmarks means that it's a not as good an approximation of the true shape overall)

#Two tests here:

  #1) Disparity of mandibles with just the landmark points (not the curves)
  #2) Disparity of mandibles with landmark points and one curve at the base of the mandible
  
#Using a reduced mandibles data file compared to my original analysis: 7 landmarks and one semilandmark curve
  #landmarks and curve are exactly the same as before (I just deleted the extra three curves)
      
library(geomorph)


source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")

#Read in the data
  setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data/mands")

#1) Data
  # Reduced data file: all landmarks and just one curve at the base of the mandible
    land <- readland.tps(file="Mands_14_03_2014_7landmarks+1bottomcurve_edited.TPS")
  
#2) Sliders
  # Reduced sliders file: just one curve at the base of the mandible
    curves <- as.matrix(read.table("Mands_14_03_2014_7landmarks+1bottomcurve_sliders_edited.NTS", header=TRUE))

#3) Taxonomy
  taxa <- read.csv("Mands_14_03_2014_Images+Specimens.csv", header=T)
#4) Specimens to remove
  rem <- read.csv("Mands_remove_spec.csv", header=T)

#--------------------------------------------
#Clean up the data
#----------------------------------------------
#Combine the taxonomic information with the array of landmark data
  combine <- list(land=land, curves=curves, ID=taxa$ID,SpecID=taxa$SpecID, Order=taxa$Order_05,
                  Fam=taxa$Family_05, Genus=taxa$Genus_05, Species=taxa$Species_05, Binom=taxa$Binomial_05)

# 1) Remove some specimens
  #1a) Remove the _sp. specimens
 sp <- which(combine$Species=="sp.")
 combine <- remove.from.list(combine, sp)
 combine <- droplevels.from.list(combine)

  #1b) remove the specimens with missing data
  matching <- matching.id(rem$SpecID, combine$SpecID)
    combine <- remove.from.list(combine, matching)
    combine <- droplevels.from.list(combine)
    
#2) Select which families to work with
  #2a) Select the tenrec and golden mole specimens only
    tc.gm <- c(which(combine$Fam=="Chrysochloridae"), which(combine$Fam=="Tenrecidae"))

    mydata <- select.from.list(combine, tc.gm)
    mydata <- droplevels.from.list(mydata)

  #2b) All of the families (just rename combine as mydata)
    #NB: remove the Notoryctidae because there's only one specimen in the family
      #mydata <- remove.from.list(combine, which(combine$Fam=="Notoryctidae"))
      #mydata <- droplevels.from.list(mydata)

#3) Option to remove the Microgale tenrecs
   mic <- which(mydata$Genus=="Microgale")

   mydata <- remove.from.list(mydata, mic)
   mydata <- droplevels.from.list(mydata)
#--------------------------------------------------
# Reduced data file with just the 7 landmarks
  
  landmarks <- array(1, dim=c(7,2,dim(mydata$land)[3]))

    for (i in 1:dim(mydata$land)[3]){
      landmarks[,,i][1:7,] <-mydata$land[,,i][c(1:7),]
    }

#Make a new data set with just these 7 landmarks and no curves

just.landmarks <- list(land=landmarks, curves=NULL, ID=mydata$ID, SpecID=mydata$SpecID, Order=mydata$Order,
                  Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

#--------------------------------------
#Procrustes superimposition
  #a) Just landmarks data
      just.landmarksGPA <- gpagen(just.landmarks$land, curves=NULL, ProcD=TRUE,)

    #List the coordinates with the taxonomic information
     Proc.co.landmarks <- list(coords=just.landmarksGPA$coords,csize=just.landmarksGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
                      Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)

  #b) Full data: 7 landmarks and one curve
      land.onecurveGPA <- gpagen(mydata$land, curves=mydata$curves, ProcD=TRUE)
      
    #List the coordinates with the taxonomic information
     Proc.co.onecurve <- list(coords=land.onecurveGPA$coords,csize=land.onecurveGPA$Csize,ID=mydata$ID,SpecID=mydata$SpecID,
                      Order=mydata$Order, Fam=mydata$Fam, Genus=mydata$Genus, Species=mydata$Species, Binom=mydata$Binom)
#-------------------------------------------------
#SPECIES AVERAGING
   #a) Just landmarks data
      #group the arrays of coordinates according to species
      just.landmarks.sps.coords <- species.coordinates(Proc.co.landmarks$coords, Proc.co.landmarks$Binom)

      #average coordinate values for each species
      just.landmarks.sps.mean <- mean.coords(just.landmarks.sps.coords)

      #list of species (same species in both data)
      binom <- just.landmarks.sps.mean$Binom
      
   #a) Full data: 7 landmarks and one curve

      onecurve.sps.coords <- species.coordinates(Proc.co.onecurve$coords, Proc.co.onecurve$Binom)
      onecurve.sps.mean <- mean.coords(onecurve.sps.coords)
      
#-----------------------------------------------
#PRINCIPAL COMPONENTS ANALYSIS
  #a) Just landmarks data
    just.landmarksPCA <- plotTangentSpace(just.landmarks.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)


  #b) Full data: 7 landmarks and 1 curve
    onecurve.PCA <- plotTangentSpace(onecurve.sps.mean$meanshape, axis1 = 1, axis2 = 2,warpgrids = TRUE, label = TRUE)

#Put the two PCA plots together
  #colour points by family
    #data frame of the unique family and binomical combinations
    sp.fam <- as.data.frame(unique(cbind(as.matrix(Proc.co.landmarks$Fam), as.matrix(Proc.co.landmarks$Binom))))
    colnames(sp.fam) <- c("Family","Binomial")

#Extract the axes for each data set
  xaxis.landmarks <- just.landmarksPCA$x[,1]
  yaxis.landmarks <- just.landmarksPCA$x[,2]

  xaxis.onecurve <- onecurve.PCA$x[,1]
  yaxis.onecurve <- onecurve.PCA$x[,2]

dev.new()
par(mfrow=c(1,2))

  plot(xaxis.landmarks,yaxis.landmarks, xlab="Landmarks: species' average PC1", ylab="Landmarks: species' average PC2",las=1,
       col=sp.fam$Family,pch=16, bty="l",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
    #draw the min,max and 0 values on the x axis
      axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
    #same for the y axis
      axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
    #add dotted lines along 0,0
      abline(0,0,h=0,v=0,lty=2,lwd=1)

  plot(xaxis.onecurve,yaxis.onecurve, xlab="Landmarks + 1curve: species' average PC1", ylab="Landmarks + 1curve: species' average PC2",las=1,
       col=sp.fam$Family,pch=16, bty="l",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
      axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
      axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
      abline(0,0,h=0,v=0,lty=2,lwd=1)

  #Seems to show higher disparity in tenrecs when the curve is included
  
#-------------------------------------------------------
#Select PC axes
  #a) Just landmarks
      PC95axes.landmarks <- selectPCaxes(just.landmarksPCA, 0.956, binom)

      #select the axes for each family
      gmolePC.landmarks <- PC95axes.landmarks[which(sp.fam$Family=="Chrysochloridae"),]
      tenrecPC.landmarks <- PC95axes.landmarks[which(sp.fam$Family=="Tenrecidae"),]
      
  #b) Landmarks + one curve
      PC95axes.onecurve <- selectPCaxes(onecurve.PCA, 0.956, binom)

      gmolePC.onecurve <- PC95axes.onecurve[which(sp.fam$Family=="Chrysochloridae"),]
      tenrecPC.onecurve <- PC95axes.onecurve[which(sp.fam$Family=="Tenrecidae"),]
      
#-------------------------------------------------------------------
#Calculate disparity
  #a) Landmarks
      #Based on PC axes

      tenrec.landmarks.sv <- PCsumvar(tenrecPC.landmarks)
      tenrec.landmarks.pv <- PCprodvar(tenrecPC.landmarks)

      tenrec.landmarks.sr <- PCsumrange(tenrecPC.landmarks)
      tenrec.landmarks.pr <- PCprodrange(tenrecPC.landmarks)

      gmole.landmarks.sv <- PCsumvar(gmolePC.landmarks)
      gmole.landmarks.pv <- PCprodvar(gmolePC.landmarks)

      gmole.landmarks.sr <- PCsumrange(gmolePC.landmarks)
      gmole.landmarks.pr <- PCprodrange(gmolePC.landmarks)

      #Based on sum of squared distances (Zeldich 2012)
        #interlandmark distance: compare each species to the overall mean shape of all species
          ild.dist.landmarks <- dist.to.ref(just.landmarks.sps.mean$meanshape)
        #create a matrix version as well with species as rownames (need it for permutation tests later in the script)
          ild.dist.landmarks.mat <- as.matrix(ild.dist.landmarks)
          rownames(ild.dist.landmarks.mat) <- binom

          tenrec.landmarks.ild <- ild.dist.landmarks[which(sp.fam$Fam == "Tenrecidae")]
          tenrecMD.landmarks<- ZelditchMD(tenrec.landmarks.ild)

          gmole.landmarks.ild <- ild.dist.landmarks[which(sp.fam$Fam == "Chrysochloridae")]
          gmoleMD.landmarks <- ZelditchMD(gmole.landmarks.ild)
          
    #Summary table

  disp.landmarks <- matrix(NA,nrow=5, ncol=2)
    rownames(disp.landmarks) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    colnames(disp.landmarks) <- c("Tenrec.obs","Gmole.obs")
    disp.landmarks[1,] <- c(tenrec.landmarks.sv, gmole.landmarks.sv)
    disp.landmarks[2,] <- c(tenrec.landmarks.pv, gmole.landmarks.pv)
    disp.landmarks[3,] <- c(tenrec.landmarks.sr, gmole.landmarks.sr)
    disp.landmarks[4,] <- c(tenrec.landmarks.pr, gmole.landmarks.pr)
    disp.landmarks[5,] <- c(tenrecMD.landmarks, gmoleMD.landmarks)
    
  #B) Landmarks and one curve
      #Based on PC axes

      tenrec.onecurve.sv <- PCsumvar(tenrecPC.onecurve)
      tenrec.onecurve.pv <- PCprodvar(tenrecPC.onecurve)

      tenrec.onecurve.sr <- PCsumrange(tenrecPC.onecurve)
      tenrec.onecurve.pr <- PCprodrange(tenrecPC.onecurve)

      gmole.onecurve.sv <- PCsumvar(gmolePC.onecurve)
      gmole.onecurve.pv <- PCprodvar(gmolePC.onecurve)

      gmole.onecurve.sr <- PCsumrange(gmolePC.onecurve)
      gmole.onecurve.pr <- PCprodrange(gmolePC.onecurve)

      #Based on sum of squared distances (Zeldich 2012)
        #interlandmark distance: compare each species to the overall mean shape of all species
          ild.dist.onecurve <- dist.to.ref(onecurve.sps.mean$meanshape)
        #create a matrix version as well with species as rownames (need it for permutation tests later in the script)
          ild.dist.onecurve.mat <- as.matrix(ild.dist.onecurve)
          rownames(ild.dist.onecurve.mat) <- binom

          tenrec.onecurve.ild <- ild.dist.onecurve[which(sp.fam$Fam == "Tenrecidae")]
          tenrecMD.onecurve<- ZelditchMD(tenrec.onecurve.ild)

          gmole.onecurve.ild <- ild.dist.onecurve[which(sp.fam$Fam == "Chrysochloridae")]
          gmoleMD.onecurve <- ZelditchMD(gmole.onecurve.ild)

    #Summary table

  disp.onecurve <- matrix(NA,nrow=5, ncol=2)
    rownames(disp.onecurve) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    colnames(disp.onecurve) <- c("Tenrec.obs","Gmole.obs")
    disp.onecurve[1,] <- c(tenrec.onecurve.sv, gmole.onecurve.sv)
    disp.onecurve[2,] <- c(tenrec.onecurve.pv, gmole.onecurve.pv)
    disp.onecurve[3,] <- c(tenrec.onecurve.sr, gmole.onecurve.sr)
    disp.onecurve[4,] <- c(tenrec.onecurve.pr, gmole.onecurve.pr)
    disp.onecurve[5,] <- c(tenrecMD.onecurve, gmoleMD.onecurve)
    
#-------------------------------------------------------------------------------
#Permutation tests for significant differences in disparity

  #a) Landmarks data
      #Pairwise observed differences in disparity among all family groups(just tenrecs vs. golden moles but could be scaled up)
        obs.diff.landmarks.sv <- group.pair.diff(group.identity=colnames(disp.landmarks), group.values=disp.landmarks[1,])
        obs.diff.landmarks.pv <- group.pair.diff(group.identity=colnames(disp.landmarks), group.values=disp.landmarks[2,])
        obs.diff.landmarks.sr <- group.pair.diff(group.identity=colnames(disp.landmarks), group.values=disp.landmarks[3,])
        obs.diff.landmarks.pr <- group.pair.diff(group.identity=colnames(disp.landmarks), group.values=disp.landmarks[4,])
        obs.diff.landmarks.md <- group.pair.diff(group.identity=colnames(disp.landmarks), group.values=disp.landmarks[5,])
  
      #Permutation tests for significant differences in disparity between tenrecs and golden moles
        #Permutation tests
          tc.gm.perm.landmarks.sv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.landmarks, PCsumvar)
          tc.gm.perm.landmarks.pv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.landmarks, PCprodvar)
          tc.gm.perm.landmarks.sr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.landmarks, PCsumrange)
          tc.gm.perm.landmarks.pr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.landmarks, PCprodrange)
          tc.gm.perm.landmarks.md <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, ild.dist.landmarks.mat, ZelditchMD)

          #test for significant differences
            tc.gm.landmarks.pvalue.sv <- pvalue.dist(tc.gm.perm.landmarks.sv, obs.diff.landmarks.sv[1,3])
            tc.gm.landmarks.pvalue.pv <- pvalue.dist(tc.gm.perm.landmarks.pv, obs.diff.landmarks.pv[1,3])
            tc.gm.landmarks.pvalue.sr <- pvalue.dist(tc.gm.perm.landmarks.sr, obs.diff.landmarks.sr[1,3])
            tc.gm.landmarks.pvalue.pr <- pvalue.dist(tc.gm.perm.landmarks.pr, obs.diff.landmarks.pr[1,3])
            tc.gm.landmarks.pvalue.md <- pvalue.dist(tc.gm.perm.landmarks.md, obs.diff.landmarks.md[1,3])

  #b) Landmarks + onecurve data
      #Pairwise observed differences in disparity among all family group
        obs.diff.onecurve.sv <- group.pair.diff(group.identity=colnames(disp.onecurve), group.values=disp.onecurve[1,])
        obs.diff.onecurve.pv <- group.pair.diff(group.identity=colnames(disp.onecurve), group.values=disp.onecurve[2,])
        obs.diff.onecurve.sr <- group.pair.diff(group.identity=colnames(disp.onecurve), group.values=disp.onecurve[3,])
        obs.diff.onecurve.pr <- group.pair.diff(group.identity=colnames(disp.onecurve), group.values=disp.onecurve[4,])
        obs.diff.onecurve.md <- group.pair.diff(group.identity=colnames(disp.onecurve), group.values=disp.onecurve[5,])

      #Permutation tests for significant differences in disparity between tenrecs and golden moles
        #Permutation tests
          tc.gm.perm.onecurve.sv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.onecurve, PCsumvar)
          tc.gm.perm.onecurve.pv <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.onecurve, PCprodvar)
          tc.gm.perm.onecurve.sr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.onecurve, PCsumrange)
          tc.gm.perm.onecurve.pr <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, PC95axes.onecurve, PCprodrange)
          tc.gm.perm.onecurve.md <- perm.diff.two.groups(1000, "Chrysochloridae", "Tenrecidae", sp.fam, ild.dist.onecurve.mat, ZelditchMD)

          #test for significant differences
            tc.gm.onecurve.pvalue.sv <- pvalue.dist(tc.gm.perm.onecurve.sv, obs.diff.onecurve.sv[1,3])
            tc.gm.onecurve.pvalue.pv <- pvalue.dist(tc.gm.perm.onecurve.pv, obs.diff.onecurve.pv[1,3])
            tc.gm.onecurve.pvalue.sr <- pvalue.dist(tc.gm.perm.onecurve.sr, obs.diff.onecurve.sr[1,3])
            tc.gm.onecurve.pvalue.pr <- pvalue.dist(tc.gm.perm.onecurve.pr, obs.diff.onecurve.pr[1,3])
            tc.gm.onecurve.pvalue.md <- pvalue.dist(tc.gm.perm.onecurve.md, obs.diff.onecurve.md[1,3])

######################################
#Add the permutation and pvalues to the summary tables

disp.landmarks.summary <- matrix(NA, nrow=5, ncol=6)
    rownames(disp.landmarks.summary) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    colnames(disp.landmarks.summary) <- c("tenrec.obs","gmole.obs", "obs.diff", "perm.min", "perm.max", "pvalue")
  disp.landmarks.summary[,1:2] <- disp.landmarks
  disp.landmarks.summary[,3] <- c(obs.diff.landmarks.sv[1,3], obs.diff.landmarks.pv[1,3], obs.diff.landmarks.sr[1,3], obs.diff.landmarks.pr[1,3], obs.diff.landmarks.md[1,3])
  disp.landmarks.summary[,4] <- c(min(tc.gm.perm.landmarks.sv), min(tc.gm.perm.landmarks.pv), min(tc.gm.perm.landmarks.sr), min(tc.gm.perm.landmarks.pr), min(tc.gm.perm.landmarks.md))
  disp.landmarks.summary[,5] <- c(max(tc.gm.perm.landmarks.sv), max(tc.gm.perm.landmarks.pv), max(tc.gm.perm.landmarks.sr), max(tc.gm.perm.landmarks.pr), max(tc.gm.perm.landmarks.md))
  disp.landmarks.summary[,6] <- c(tc.gm.landmarks.pvalue.sv,  tc.gm.landmarks.pvalue.pv,  tc.gm.landmarks.pvalue.sr,  tc.gm.landmarks.pvalue.pr,  tc.gm.landmarks.pvalue.md)
  
disp.onecurve.summary <- matrix(NA, nrow=5, ncol=6)
    rownames(disp.onecurve.summary) <- c("SumVar","ProdVar","SumRange","ProdRange", "ZelditchMD")
    colnames(disp.onecurve.summary) <- c("tenrec.obs","gmole.obs", "obs.diff", "perm.min", "perm.max", "pvalue")
  disp.onecurve.summary[,1:2] <- disp.landmarks
  disp.onecurve.summary[,3] <- c(obs.diff.onecurve.sv[1,3], obs.diff.onecurve.pv[1,3], obs.diff.onecurve.sr[1,3], obs.diff.onecurve.pr[1,3], obs.diff.onecurve.md[1,3])
  disp.onecurve.summary[,4] <- c(min(tc.gm.perm.onecurve.sv), min(tc.gm.perm.onecurve.pv), min(tc.gm.perm.onecurve.sr), min(tc.gm.perm.onecurve.pr), min(tc.gm.perm.onecurve.md))
  disp.onecurve.summary[,5] <- c(max(tc.gm.perm.onecurve.sv), max(tc.gm.perm.onecurve.pv), max(tc.gm.perm.onecurve.sr), max(tc.gm.perm.onecurve.pr), max(tc.gm.perm.onecurve.md))
  disp.onecurve.summary[,6] <- c(tc.gm.onecurve.pvalue.sv,  tc.gm.onecurve.pvalue.pv,  tc.gm.onecurve.pvalue.sr,  tc.gm.onecurve.pvalue.pr,  tc.gm.onecurve.pvalue.md)



#Significant difference in the ZelditchMD for both data sets
  #Otherwise there's no significant difference in tenrec and golden mole disparity
  #So the differences in my full analysis come from more variable muscle attachment shapes in golden moles
    #Now the question is why!
#---------------------------------------------    
#Output
    setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/output/shape_data/mands")

#PCA plot for the mandibles data with one curve
  #All tenrecs and golden moles
  #pdf(file="Mandibles_trc+gmole_onecurve_PCA.pdf")
  
  #Non-Microgale tenrecs and golden moles  
  pdf(file="Mandibles_nonmic_trc+gmole_onecurve_PCA.pdf")  
    
    plot(xaxis.onecurve,yaxis.onecurve, xlab="Species' average PC1", ylab="Species' average PC2",las=1,
       col=sp.fam$Family,pch=16, bty="l",cex.lab=1,cex=1.2, xaxt="n",yaxt="n")
      axis(side=1,at=c(-0.1,0,0.1),las=1,cex=1.3)
      axis(side=2,at=c(-0.06,0,0.08),las=1,cex=1.3)
      abline(0,0,h=0,v=0,lty=2,lwd=1)
  dev.off()
  
#Summary table of disparity calculations
  #All terecs and golden moles
  #write.table(file="Mandibles_trc+gmole_onecurve_disp_summary.txt",disp.onecurve.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)
  
  #Non-Microgale tenrecs and golden moles
  write.table(file="Mandibles_nonmic_trc+gmole_onecurve_disp_summary.txt",disp.onecurve.summary,col.names=T, row.names=T,sep="\t",quote=F,append=FALSE)  

#generate a latex table directly: then copy and paste it into a text file
library(xtable)
xtable(disp.onecurve.summary)