#Functions used in my analysis of morphological diversity in tenrecs and golden moles
  #Code came from my original functions scripts but it's neater to just have the functions that I'm actually using
  #Revised analysis based on comparing mean distances to centroid in each family
  
#1) General functions
    # matching.id
    # remove.from.list
    # select.from.list
    # droplevels.from.list
    # anova.frp
  
#2) Shape data
    # species.coordinates
    # mean.coords
    # selectPCaxes
    
#3) Comparing morphological diversity in two groups
    # euc.dist.cent
    # group.diff
    # pvalue.dist
    
####################################
#1) General functions
######################
#Function to find the ID (row) numbers that match a subset of specimens
  #subset.col and main.col are the relevant columns in each data set
  matching.id <- function (subset.col, main.col){
    id.list <- rep(NA, length(subset.col))

      for (i in 1:length (subset.col)){
        id.list[i] <- grep(subset.col[i], main.col)
      }
    return(id.list)
  }
#------------------------------------------------
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

#---------------------------------
#Extract the f value, r squared and p value from an anova table
  anova.frp <- function(anova.object){
    anova.frp <- matrix(NA, nrow=1, ncol=3)
    colnames(anova.frp) <- c("F.Model", "R2", "p")
      anova.frp[1,1] <- anova.object$aov.tab$F.Model[1]
      anova.frp[1,2] <- anova.object$aov.tab$R2[1]
      anova.frp[1,3] <- anova.object$aov.tab$Pr[1]
    return(anova.frp)
  }
##########################################
#2) Shape data
##########################################
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
#---------------------------------------------------
#Function to select specific PC axes from a pcaresults object
  #select based on a threshold amount of cumulative variation explained
    #avoids selecting just single axes
    #just selects the axes that meet the threshold (original version added 1 extra axis on top of that)

  selectPCaxes <- function(pcaresults, threshold, species){
    if (pcaresults$importance[3,1] > threshold){
        no.of.axes <- 1
    } else {
        no.of.axes <- length(which(pcaresults$importance[3,] <= threshold))
        }
        PCaxes <- pcaresults$x[,1:(no.of.axes)]
        rownames(PCaxes) <- species
    return(PCaxes)
  }



##########################################
#3) Comparing the diversity of two groups
##########################################

#Function to calculate mean euclidean distances from the centroid
  euc.dist.cent <- function(PCdata){
    #Centroid (mean score of each PC axis)
    centroid <- NULL
      for(i in 1:ncol(PCdata)){
        centroid[i] <- mean(PCdata[,i])
      }

    #Euclidean distances to the centroid
    cent.dist <- NULL
      for (j in 1:nrow(PCdata)){
        cent.dist[j] <- dist(rbind(PCdata[j,], centroid), method="euclidean")
      }
    return(cent.dist)
  }


#---------------------------------------------------
#Function for permutation tests of significant differences between two groups
  #Modified from Steve Wang's email on 10/06/2014

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

#-------------------------
#Function to calculate a p value from a distribution
  #distribution is the set of re-sampled measurements
  #obs.val is the observed value of the test statistic

pvalue.dist <- function(distribution,obs.val){
  low.diff <- distribution[which(distribution <= obs.val)]
    length(low.diff)
      lowerp <- (length(low.diff))/(length(distribution))
        higherp <- (1-lowerp)

          if (lowerp < higherp){
            p <- lowerp

          } else {
          p <- higherp
        }
    #print both the lower and higher p values
    cat ("lowerp",lowerp,"\n")
    cat ("higherp",higherp,"\n")
    #but return the actual p value
  return(p)
}