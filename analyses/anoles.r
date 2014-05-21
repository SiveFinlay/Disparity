#Test disparity code with anolis data
#Data from dryad repositories of Harmon et al 2010 and Thomas et al 2009

library(vegan)

setwd("C:/Users/sfinlay/Desktop/Thesis/Disparity/data")

source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/Disparity_general_functions.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/DisparityFunctions_Variance_Range.r")
source("C:/Users/sfinlay/Desktop/Thesis/Disparity/functions/PvalueFunction_FromDistribution.r")
########################################################
#Extract and combine common species from each data set

#Anole data from Harmon et al 2010
  anole.harmon <- read.table("anolisraw.dat", header=T, sep = "\t", na.strings=".")

#Anole data from Thomas et al 2009
  anole.thomas <- read.table("Thomas2009_anolisdata.txt", header=TRUE, sep="\t")
  anole.thomas.species <- split.binom(as.character(anole.thomas$Species))

#list of species in each data
  an.th <- unique(anole.thomas.species$Species)
  an.har <- unique(as.character(anole.harmon$NAME))

#find the common species in the two data sets
  common <- common.character(an.th, an.har)

#select those common species from each data set
  an.th.com.id <- NULL
    for (i in 1:length(common)){
      an.th.com.id[i] <- which (common[i] == anole.thomas.species$Species)
    }

  an.th.com <- anole.thomas[an.th.com.id,] 
  #pick out the relevant columns
  anole.th <- cbind(anole.thomas.species[an.th.com.id,], island=an.th.com$Island_type, 
                    ecomorph=an.th.com$ecomorph, femaleSVL=an.th.com$Female_SVL, maleSVL=an.th.com$Male_SVL)
    
  an.har.com.id <- NULL
    for (i in 1:length(common)){
      an.har.com.id[i] <- which (common[i] == anole.harmon$NAME)
    }
  
  an.har.com <- anole.harmon[an.har.com.id,]

#Combine into one data frame
  anole.data <- droplevels(cbind(anole.th[,3:6], an.har.com[,2:7]))  
  rownames(anole.data) <- anole.th$Species
#-------------------------------------------------------------------------
#Compare disparity of species from large and small islands

#PCA of all data
anole.PCA <- prcomp(anole.data[,3:10], data=anole.PCA, scale=TRUE )

#Select PC axes that account for 95% of the variation
  anolePC95 <- selectPCaxes.prcomp(anole.PCA, 0.956)

#Subset the data into large and small island groups
  large <- anolePC95[which(anole.data$island == "LargeIsland"),]
  small <- anolePC95[which(anole.data$island == "SmallIsland"),]
  
#Calculate disparity
  #Based on PC axes
  large.sv <- PCsumvar(large)
  large.pv <- PCprodvar(large)
  large.sr <- PCsumrange(large)
  large.pr <- PCprodrange(large)
  
  small.sv <- PCsumvar(small)
  small.pv <- PCprodvar(small)
  small.sr <- PCsumrange(small)
  small.pr <- PCprodrange(small)

#Put the disparity calculations into a single table
  disp.island <- matrix(NA,nrow=2, ncol=4)
    rownames(disp.island) <- c("Large","Small")
    colnames(disp.island) <- c("SumVar","ProdVar","SumRange","ProdRange")
    disp.island[,1] <- c(large.sv,small.sv)
    disp.island[,2] <- c(large.pv,small.pv)
    disp.island[,3] <- c(large.sr,small.sr)
    disp.island[,4] <- c(large.pr,small.pr)

#Compare the two groups
  #1) Distance matrix

  #Euclidean distance matrix of all of the species
    Euc.dist <- as.matrix(dist(anolePC95, method="euclidean", diag=FALSE,upper=FALSE))

#npMANOVA of the distance matrix separated by island group
  dist.man.isl <- adonis(Euc.dist~anole.data$island , data=anole.data, permutations=9999,method="euclidean")
    #extract the f, r2 and p values
      dist.man.isl.frp <- anova.frp(dist.man.isl)    

  #2) PC axes
  #NPMANOVA of the PC axes  
  PC.man.isl <- adonis(anolePC95~anole.data$island , data=anole.data, permutations=999, method="euclidean")
    #extract the f, r2 and p values
      PC.man.isl.frp <- anova.frp(PC.man.isl)

#Compare the distance and PC npMANOVA results    #no significant difference
  manova.res.isl <- rbind(dist.man.isl.frp, PC.man.isl.frp)
  rownames(manova.res.isl) <-c ("dist.man", "PC.man")
    #NB: Results from Thomas et al 2009: rates of diversification on small islands are higher than large islands
        # i.e. they didn't compare current patterns of disparity explicitly

#------------------------------------------------------------------------------
#Same approach but compare disparity of ecomorph groups  
  trunkground <- anolePC95[which(anole.data$ecomorph == "TrunkGround"),]
  trunkcrown <- anolePC95[which(anole.data$ecomorph == "TrunkCrown"),]
  #combine the two
  trunk <- rbind(trunkground, trunkcrown)
  #combined classification
  trunk.class <- as.data.frame(cbind(species = rownames(trunk), 
                   ecomorph = c(rep("ground", nrow(trunkground)), rep("crown", nrow(trunkcrown)))))
  
#Calculate disparity
  #Based on PC axes
  trunkground.sv <- PCsumvar(trunkground)
  trunkground.pv <- PCprodvar(trunkground)
  trunkground.sr <- PCsumrange(trunkground)
  trunkground.pr <- PCprodrange(trunkground)
  
  trunkcrown.sv <- PCsumvar(trunkcrown)
  trunkcrown.pv <- PCprodvar(trunkcrown)
  trunkcrown.sr <- PCsumrange(trunkcrown)
  trunkcrown.pr <- PCprodrange(trunkcrown)

#Put the disparity calculations into a single table
  disp.eco <- matrix(NA,nrow=2, ncol=4)
    rownames(disp.eco) <- c("TrunkGround","TrunkCrown")
    colnames(disp.eco) <- c("SumVar","ProdVar","SumRange","ProdRange")
    disp.eco[,1] <- c(trunkground.sv,trunkcrown.sv)
    disp.eco[,2] <- c(trunkground.pv,trunkcrown.pv)
    disp.eco[,3] <- c(trunkground.sr,trunkcrown.sr)
    disp.eco[,4] <- c(trunkground.pr,trunkcrown.pr)

#Compare the two groups

  #Euclidean distance matrix of all of the species
    Euc.dist.trunk <- as.matrix(dist(trunk, method="euclidean", diag=FALSE,upper=FALSE))

#npMANOVA of the distance matrix separated by island group
  dist.man.eco <- adonis(Euc.dist.trunk~trunk.class$ecomorph , data=trunk.class, permutations=9999,method="euclidean")
    #extract the f, r2 and p values
      dist.man.eco.frp <- anova.frp(dist.man.eco)    

  #NPMANOVA of the PC axes  
  PC.man.eco <- adonis(trunk~trunk.class$ecomorph , data=trunk.class, permutations=999, method="euclidean")
    #extract the f, r2 and p values
      PC.man.eco.frp <- anova.frp(PC.man.eco)

#Compare the distance and PC npMANOVA results    #no significant difference
  manova.res.eco <- rbind(dist.man.eco.frp, PC.man.eco.frp)
  rownames(manova.res.eco) <-c ("dist.man", "PC.man")      
      #significant difference in disparity between trunkground and trunkcrown ecomorphs
