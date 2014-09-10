#Start date: 06/05/2014


#Useful functions from my re-writing of the disparity analysis
  #Use this as a source file for all of my analyses, not just the disparity ones
  #Includes functions that are generally useful for data manipulation and shape data

#1) General functions
    #matching.id 
    #subset.matrix
    #anova.frp
    #common.character
    #group.pair.diff
    #group.diff (cf Steve Wang)
    #perm.diff.two.groups (wrapper function for group.diff)
    
#2) Dealing with lists
      #remove.from.list
      #select.from.list
      #droplevels.from.list
      #calc.each.list
      #calc.each.array
      #list.arrays.to.matrices
      #sorted.list
      
#3) Dealing with shape data
      #species.coordinates
      #mean.coordinates
      #selectPCaxes
      #selectPCaxes.prcomp
      #Proc.dist.within

#4) Resampling (rarefaction)
      #resample.data
      #resample.data.replace
      #mean.fun
      #boot.mean.1000
      #boot.95.min.confidence
      #boot.95.max.confidence

      
#5) Dealing with phylogenies
      #tree.only
      #remove.missing.species.tree
        #(remove.missing.species.phy and remove.missing.species.multiphy)
      #split.binom
      #add.species.to.MRCA

#6) Plotting functions
      #arrow.to.x.point

#******************************************
#1) General functions
#******************************************
#Function to find the ID (row) numbers that match a subset of specimens
  #subset.col and main.col are the relevant columns in each data set
  matching.id <- function (subset.col, main.col){
    id.list <- rep(NA, length(subset.col))
 
      for (i in 1:length (subset.col)){
        id.list[i] <- grep(subset.col[i], main.col)
      }
    return(id.list)
  }
  
#--------------------------------------
#Function to select particular rows from a matrix
  subset.matrix <- function(mydata, subset.col, criteria){
    set <- mydata[which(subset.col == criteria),]
  }
  
#------------------------------------
#Extract the f value, r squared and p value from an anova table
  anova.frp <- function(anova.object){
    anova.frp <- matrix(NA, nrow=1, ncol=3)
    colnames(anova.frp) <- c("F.Model", "R2", "p")
      anova.frp[1,1] <- anova.object$aov.tab$F.Model[1]
      anova.frp[1,2] <- anova.object$aov.tab$R2[1]
      anova.frp[1,3] <- anova.object$aov.tab$Pr[1]
    return(anova.frp)
  }
  
#--------------------------------------------------------  
#Function to find common characters in two lists
  common.character <- function (long.list, short.list){
    tmp<-NULL
    value<-NULL
    common<-NULL
      for (i in 1: length(long.list)){
        tmp<-(long.list[i] == short.list)
        value<-which(tmp == TRUE)
        common<-c(common, short.list[value])
      }
      return(common)
  }
#---------------------------------------------------------
#Function to calculate pairwise differences among single values for different groups
  #Returns a data frame listing the two groups and the differences in their values

  group.pair.diff <- function (group.identity, group.values){
    group.pair <- combn(group.identity,2)
    value.pair.diff <- diff(combn(group.values, 2))
    pair.diff <- matrix(NA, nrow=ncol(group.pair), ncol=3)
      colnames(pair.diff) <- c("group1", "group2", "group2-group1") 
        for (i in 1:ncol(group.pair)){
          pair.diff[i,1] <- group.pair[,i][1]
          pair.diff[i,2] <- group.pair[,i][2]
          pair.diff[i,3] <- value.pair.diff[i]
        }  
    pair.diff <- as.data.frame(pair.diff)
    pair.diff[,3] <- as.numeric(as.character(pair.diff[,3]))
    return(pair.diff)
  }
  
#---------------------------------------------------
#Function for permutation tests of significant differences between two groups
  #Modified from Steve Wang's email on 10/06/2014
  #Modified again on 16/07/2014 to make the function applicable to more than two groups

  
  group.diff <- function (numreps, two.groups, mydata, test.statistic){
    results <- rep(NA,numreps)
    for (rep in 1:numreps) {
      shufgroup <- sample(two.groups)
        shuff.group1 <- mydata[which(shufgroup==(levels(two.groups)[1])),]
        shuff.group2 <- mydata[which(shufgroup==(levels(two.groups)[2])),]
          shuffdiff <- test.statistic(shuff.group1) - test.statistic(shuff.group2)
      results[rep] <- shuffdiff
    }
    return(results)
  } 

#Wrapper function for using group.diff; select two family groups and corresponding data from larger data objects

perm.diff.two.groups <- function (numreps, fam1, fam2, sp.fam.data, mydata, test.statistic){
  #list of species names
  mygroups <- droplevels(sp.fam.data[c(which(sp.fam.data$Family == fam1), (which(sp.fam.data == fam2))),])
  #select those species from mydata
  mygroups.id <- NULL
    for (i in 1:length(mygroups$Binomial)){
      mygroups.id[i] <- which(rownames(mydata) == mygroups$Binomial[i])
    }
  mygroups.data <- mydata[mygroups.id,,drop=FALSE]    #drop=FALSE preserves the matrix class of the selected object
  #permutation test for significant difference in a test statistic between the two groups
  perm.group <- group.diff(numreps, mygroups$Family, mygroups.data, test.statistic)
}


#****************************************
#2) Dealing with lists
#****************************************
#Function to remove a vector of numbers from non matrix objects in a list
  remove.from.list <- function(mylist, remove.id){
    newlist <- as.list(rep(NA,length(mylist)))
    names(newlist) <- names(mylist)

      for (i in 1: length(mylist)){
        if(class(mylist[[i]]) == "matrix"){
          newlist[[i]]<-mylist[[i]]  #don't remove anything from matrix elements (curves in my landmark data)

        } else {        
          if(class(mylist[[i]]) == "array"){
            newlist[[i]] <- mylist[[i]][,,-remove.id] #remove from the third dimension of an array (landmark coordinates)

          } else {  
            newlist[[i]] <- mylist[[i]][-remove.id] #remove from all other elements in the list
            }
          }
      }
      return(newlist)
  }

#---------------------------------------------
#Function to select specific ID numbers from non-matrix objects in a list  
  select.from.list <- function(mylist, select.id){
    newlist <- as.list(rep(NA,length(mylist))) 
    names(newlist) <- names(mylist)       
      
      for (i in 1: length(mylist)){
       if (class(mylist[[i]]) == "matrix"){
         newlist[[i]] <- mylist[[i]] #leave matrix elements in the list unchanged (curves in the landmark data) 
      
       } else {         
         if(class(mylist[[i]]) == "array"){
           newlist[[i]] <- mylist[[i]][,,select.id] #select from the third dimension of an array (landmark coordinates)
      
         } else { 
           newlist[[i]] <- mylist[[i]][select.id]   #select from all other elements in the list
           }
         }
      }
      return(newlist)
  }
  
#------------------------------------------------------------------
#Function to drop unused levels from list elements that are factors and characters
  droplevels.from.list <- function(mylist){
    newlist <- as.list(rep(NA,length(mylist)))
    names(newlist) <- names(mylist)        
      
      for (i in 1: length(mylist)) {
        if (class(mylist[[i]]) != "integer" & class(mylist[[i]]) != "array" 
            & class(mylist[[i]]) != "numeric" & class(mylist[[i]]) != "matrix"){
          newlist[[i]] <- droplevels(mylist[[i]])
      
        } else {
          newlist[[i]] <- mylist[[i]]
          }
        }
        return(newlist)
  }
  
#------------------------------------------------------------------------------------
#Function to apply a calculation to each element in a list 
  calc.each.list <- function(mylist, calculation){
    new.list <- NULL

      for (i in 1:length(mylist)){
        new.list[[i]] <- calculation(mylist[[i]])
      }
     class(new.list) <- class(mylist)
     return(new.list)
  }

#----------------------------------------------------------
#Function to apply a calculation to each array within a list of arrays
  calc.each.array <- function(array.list, calculation){
    new.array.list <- NULL
 
      for (i in 1:length(array.list)){
        new.array.list[[i]] <- rep(NA, dim(array.list[[i]])[3])
      }
      #fill the empty list with a calculation for each array

        for (j in 1:length(array.list)){
          for(k in 1:dim(array.list[[j]])[3]){
            new.array.list[[j]][k] <- calculation(array.list[[j]][,,k])
          }
        }
        return(new.array.list)
  }
  
#-----------------------------------------------------------------------------------
#Function to break a list of arrays into a list of matrices    
  list.arrays.to.matrices <- function(mylist){
    new.list <- NULL

      for (i in 1:length(mylist)){
        for(m in 1:(dim(mylist[[i]])[3])){ #gives the number of arrays
          new.list[[m+(i*(dim(mylist[[i]])[3])-(dim(mylist[[i]])[3]))]] <- mylist[[i]][,,m]
                   #For 10 arrays ((dim(mylist[[i]])[3]) =10)
                    #when i is 1, m will take the values 1:10
                    #when i is 2, m will take the values 11:20
        }
      }
      return(new.list)
  }

#------------------------------------------------------------------
#Function to sort the elements in a list (from smallest to largest)
  sorted.list <- function (mylist){
    sorted.list <- NULL
      for (i in 1:length(mylist)){
        sorted.list[[i]] <- sort(mylist[[i]])
      }
    return(sorted.list)
  }  
    
#***************************************************
#3) Dealing with shape data
#***************************************************
#Function to select the Procrustes coordinates of each species
  #Inputs are an array object of the coordinates and the binomial species values that correspond to that array
  species.coordinates <- function(coords, coords.binom){
    binom.list <- unique(coords.binom)      
      #list of ID numbers for each species
      species <- as.list(rep(NA, length(binom.list)))

        for (i in 1:length(binom.list)){
          species[[i]] <- which(coords.binom == binom.list[i])
        }
      #coordinates of those ID numbers
      sps.proc.co <- as.list(rep(NA, length(binom.list)))
      names(sps.proc.co)<-binom.list

        for (j in 1:length(species)){
          sps.proc.co[[j]] <- coords[,,species[[j]]]
        }
        return (sps.proc.co)
  }

#-----------------------------------------------------------
#Function to find the average coordinates of each species
  mean.coords <- function(sps.coords){
    sps.coords.mean<-array(data=NA, dim=c(dim(sps.coords[[1]])[1],2,length(sps.coords)))
      #Select each species

        for (k in 1:length(sps.coords)){
          #get the meanshape of the aligned coordinates of that species
          for(m in 1:length(dim(sps.coords[[k]]))){
            #Only calculate the mean shape of coordinates when there's more than one set of landmarks
             if (length(dim(sps.coords[[k]])) != 2){
               sps.coords.mean[,,k] <- mshape(sps.coords[[k]])
             
             } else {
               sps.coords.mean[,,k] <- sps.coords[[k]]
               }     
           }   
        }
        sps.mean <- list(meanshape=sps.coords.mean, ID=1:length(names(sps.coords)), Binom=as.factor(names(sps.coords)))
        return (sps.mean)
  }

#-------------------------------------
#Function to select specific PC axes from a pcaresults object
  #select based on a threshold amount of cumulative variation explained and then add 1 extra axis
    #avoids selecting just single axes
      #03/06/2014 Changed function for twofamily_disparity script
              #took out $pc.summary before $importance
              #changed $pc.scores to $x
      #23/07/2014: Changed the function so that it will still select axes even if PC1 is greater than the threshold
  
  selectPCaxes <- function(pcaresults, threshold, species){                
    if (pcaresults$importance[3,1] > threshold){
        no.of.axes <- 1
    } else {
        no.of.axes <- length(which(pcaresults$importance[3,] <= threshold))      
        }
        PCaxes <- pcaresults$x[,1:(no.of.axes + 1)]   
        rownames(PCaxes) <- species
    return(PCaxes)
  }
    
#Functionto select specific PC axes from a prcomp object
  #selects based on a threshold amount of cumulative variation explained and then add 1 extra axis
  selectPCaxes.prcomp <- function (prcomp.object, threshold){
    no.of.axes <- NULL
      if(summary(prcomp.object)$importance[3,1] > threshold){
        no.of.axes <- 1
      } else {
        no.of.axes <- length(which(summary(prcomp.object)$importance[3,] <= threshold))
        }
        PCaxes <- prcomp.object$x[,1:(no.of.axes +1)]
        rownames(PCaxes) <- rownames(prcomp.object$x) 
     return(PCaxes)
   } 

#-------------------------------------------------
#Function to find the min, max and mean Procrustes distances within a set of specimens
  #Based on a Euclidean distance matrix
  #I used it for morphometrics error checking so there are multiple measures of single specimens
    #NB: The function selects the minimum non 0 number
  #(there will always be 0s in the distance matrix because there's no difference between an image and itself)   
  
  Proc.dist.within <- function (dist.mat){
    specimen <- NULL
    for (i in 1:length(unique(rownames(dist.mat)))){
      specimen[[i]] <-   which(rownames(dist.mat) == (unique(rownames(dist.mat))[i]))
    }
    specimen.dist <- NULL
    for (j in 1:length(specimen)){
      specimen.dist[[j]] <- dist.mat[specimen[[j]], specimen[[j]]]
    }
    
    dist.summary <- matrix(nrow=(length(specimen)), ncol=3)
      rownames(dist.summary) <- unique(rownames(dist.mat))
      colnames(dist.summary) <- c("non_0_min", "max", "mean")

    for (k in 1:length(specimen.dist)){
         #dist.summary[j,] <- c(range (dist.mat[specimen[[j]],specimen[[j]]]), mean(dist.mat[specimen[[j]],specimen[[j]]]))
         dist.summary[k,1] <- min(specimen.dist[[k]][which(specimen.dist[[k]] > 0)])
         dist.summary[k,2] <- max(specimen.dist[[k]])
         dist.summary[k,3] <- mean(specimen.dist[[k]])
         }
  
    return(dist.summary)
    }

#**********************************************
#4) Resampling (rarefaction)
#*********************************************
#Function to resample data without replacement for rarefaction
  resample.data <- function(mydata, samp.min, samp.max, no.replicates, no.col){
    #make a list of empty arrays first
    resample <- as.list(rep(NA, samp.max))
   
      for (i in samp.min:samp.max){
        resample[[i]] <- array(NA, dim=(c(i,no.col, no.replicates)))
      }
      #fill the empty arrays with resampled data
   
      for (j in samp.min:samp.max){
        for (k in 1:no.replicates){
          resample[[j]][,,k] <- mydata[sample(nrow(mydata), size=j, replace=FALSE),]
        }
      }
     #if samp.min is 2 then the first value of resample is NULL (didn't select multiple replicates of a single species)
      #remove that null value
      if (is.na(resample[[1]]) == TRUE){
        resample[[1]] <- NULL
      }
      return(resample)
  }

#Function to resample data with replacement  
  resample.data.replace <- function(mydata, samp.min, samp.max, no.replicates, no.col){
    #make a list of empty arrays first
    resample <- as.list(rep(NA, samp.max))
   
      for (i in samp.min:samp.max){
        resample[[i]] <- array(NA, dim=(c(i,no.col, no.replicates)))
      }
      #fill the empty arrays with resampled data
   
      for (j in samp.min:samp.max){
        for (k in 1:no.replicates){
          resample[[j]][,,k] <- mydata[sample(nrow(mydata), size=j, replace=TRUE),]
        }
      }
     #if samp.min is 2 then the first value of resample is NULL (didn't select multiple replicates of a single species)
      #remove that null value
      if (is.na(resample[[1]]) == TRUE){
        resample[[1]] <- NULL
      }
      return(resample)
  }
  
#--------------------------------------------
#Function to calculate the mean
  #non-parametric bootstrapping requires functions with at least 2 arguments
  #code for a mean function here https://stat.ethz.ch/pipermail/r-help/2008-April/160254.html
  mean.fun <- function(dat, idx) mean(dat[idx], na.rm = TRUE)

#Function to calculate mean for 1000 bootstrapped values
  boot.mean.1000 <- function (sample.values){
    boot.mean.values <- boot(data=sample.values, statistic=mean.fun, R=1000, sim="ordinary")
  }

#------------------------------------------------
#Function to return the minimum 95% confidence value from bootstrapped data
  boot.95.min.confidence <- function (boot.data){
    conf.intervals <- boot.ci(boot.data, type="norm")
    min.confidence <- conf.intervals$normal[2] 
  return (min.confidence)  
  }
  
#--------------------------------------------
#Function to return the maximum 95% confidence value from bootstrapped data
  boot.95.max.confidence <- function (boot.data){
    conf.intervals <- boot.ci(boot.data, type="norm")
    max.confidence <- conf.intervals$normal[3]
  return (max.confidence)
  }

#***************************************
#5) Dealing with phylogenies
#***************************************
#Functions to identify taxa that are in the trees but not the data

#Identify tree-only taxa in a single tree
  tree.only.phy <- function(phy,data.species){
    taxa.tree.only <- setdiff(phy$tip.label, data.species)
    }

#Identify tree-only taxa in multiple trees
  tree.only.multiphy <- function(multiphy,data.species){
    taxa.tree.only <- NULL
      for (i in 1:length(multiphy)){
        taxa.tree.only[[i]] <- setdiff(multiphy[[i]]$tip.label, data.species)
          }
          return(taxa.tree.only)
  }
    
#Wrapper function to select tree-only taxa from either a phy or multiphy object    
  tree.only <- function(phy,data.species){
    taxa.tree.only <- NULL
      if(class(phy) == "phylo"){
        taxa.tree.only <- tree.only.phy(phy, data.species)
      } else {  
        taxa.tree.only <- tree.only.multiphy(phy, data.species)
        }
        return(taxa.tree.only)
  }
    
#---------------------------------------------    
#Functions to remove missing species from trees

#Remove missing species from a single tree
  remove.missing.species.phy <- function(phy, missing.species) {
    new.trees <- drop.tip(phy, missing.species)
  }

#Remove missing species from multiple trees
  remove.missing.species.multiphy <- function(multiphy, missing.species){
    new.trees<-NULL
      for (i in 1:length(multiphy)){
        new.trees[[i]] <- drop.tip(multiphy[[i]], missing.species[[i]])
      }
      return(new.trees)
  }

#Wrapper function to remove missing species from either a phy or MultiPhy object
  remove.missing.species.tree <- function(mytrees, missing.species){
     new.trees <- NULL
       if (class(mytrees) == "phylo"){
         new.trees <- remove.missing.species.phy(mytrees, missing.species)
 
       } else {
         new.trees <- remove.missing.species.multiphy(mytrees, missing.species)
         }
         return(new.trees)
  }
  
#----------------------------------------------
#Function to split Binomial names into a data frame of genus and species names
  split.binom <- function(binom.list){
    sps.split <- strsplit(binom.list, split="_")  
    split.names <- as.data.frame(matrix(NA, nrow=length(sps.split), ncol=2))
      colnames(split.names) <- c("Genus", "Species")
        
        for (i in 1:nrow(split.names)){
          split.names[i,1] <- sps.split[[i]][1]
            split.names[i,2] <- sps.split[[i]][2]
        }
        return(split.names)
    }

#-------------------------------------------------------
#Function to add a species at random to the most recent common ancestor of a list of species
  add.species.to.MRCA <- function(mytrees, species.in.tree, species.to.add){
  #find the most recent common ancestor of the species
    anc<-NA
      for (i in 1:length(mytrees)){
        anc[i] <- findMRCA(tree=mytrees[[i]], tips=species.in.tree, type="node")
      }
  #bind the new species to that ancestral node
    newtrees<-mytrees
      for (j in 1:length(newtrees)){
        newtrees[[j]] <- bind.tip(tree=mytrees[[j]], tip.label=species.to.add, edge.length=NULL, where=anc[j], position=0)
      }
      return(newtrees)
  }
                      
###########################################################
#6) PLOTTING FUNCTIONS
#############################################################

#Function to add an arrow to a histogram which points to a particular x value
  arrow.to.x.point <- function (histog, x.value, fraction.of.yaxis, line.fraction.of.yaxis, height.above.xaxis, head.length, colour, line.width) {
    point <- c(x.value, (max(histog$counts)/fraction.of.yaxis))
    line.length <- max(histog$counts)/line.fraction.of.yaxis
      
      arrows((point[1]), (point[2] + line.length), (point[1]), (point[2] + height.above.xaxis), 
            length=head.length, col=colour, lwd=line.width, code=2) 
  }
