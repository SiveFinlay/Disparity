#06/05/2014
#07/05/2014
#12/05/2014

#Useful functions from my re-writing of the disparity analysis

#1) General functions
    #matching.id 
    #subset.matrix
    #anova.frp
    
#2) Dealing with lists
      #remove.from.list
      #select.from.list
      #droplevels.from.list
      #calc.each.list
      #calc.each.array
      #list.matrix.to.array
      
#3) Dealing with shape data
      #species.coordinates
      #mean.coordinates
      #selectPCaxes

#4) Resampling (rarefaction)
      #resample.data
      
#5) Dealing with phylogenies
      #tree.only
      #remove.missing.species.tree
        #(remove.missing.species.phy and remove.missing.species.multiphy)
      #split.binom
      #add.species.to.MRCA

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
#Function to turn a list of matrices into a list of arrays    
  list.matrix.to.array <- function(mylist){
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
  #select based on a threshold and then add 1 extra axis
    #avoids selecting just single axes
  selectPCaxes <- function(pcaresults, threshold, species){
    no.of.axes <- length(which(pcaresults$pc.summary$importance[3,] <= threshold))
      PCaxes <- pcaresults$pc.scores[,1:(no.of.axes + 1)]
      rownames(PCaxes) <- species
    return(PCaxes)
  }
    
#**********************************************
#4) Resampling (rarefaction)
#*********************************************
#Function to resample data for rarefaction
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
                      