#13/05/2014 
  #Updated the code following Google's R style guide
  #Shortened functions
    #Original code for disparity functions in the DisparityPractice_16_01_14_using_SkLat_08_11_13 script

#Functions
  #1) Based on PC axes
    #PCvariance
    #PCsumvar
    #PCprodvar
  
    #PCrange
    #PCsumrange
    #PCprodrange

  #2) Based on interlandmark distances
    #ild
    #dist.to.ref
    #ZelditchMD

################################################ 
#1) Based on PC axes
################################################
#Calculate the variance of each PC axis
  PCvariance <- function(PCdata){
    PCmatrix <- as.matrix(PCdata)
    PCvar <- NULL
      for (i in 1:ncol(PCmatrix)){
        PCvar[i] <- var(PCmatrix[,i])
      }
      return(PCvar)
  }

#------------------------------------------------------------
#Caculate the sum of variance 
  PCsumvar <- function(PCdata){
   PCvar <- PCvariance(PCdata)
    sumvar <- NULL
      if (length(PCvar) > 1){
        sumvar <- sum (PCvar)
      } else {
        sumvar <- 0 #single axes can't have a sum of variance
        }
    return(sumvar)
  }

#------------------------------------------------------------
#Calculate the product of variance
  PCprodvar <- function(PCdata){
    PCvar <- PCvariance(PCdata)
      prodvar<-NULL
        if (length(PCvar) > 1) {
            prodvar <-prod(PCvar)
        } else {
          prodvar <- 0     #single axes can't have a product of variance
          } 
          #divide the product by the root of the number of axes used to reduce the dimensionality of the answer (cf Wills + Brusatte)
          prodvar.scaled<-prodvar^(1/(length(PCvar)))
          return(prodvar.scaled)
  }

#------------------------------------------------------------
#Calculate the range of each PC axis
  PCrange <- function(PCdata){
    PCmatrix <- as.matrix(PCdata)
      PCrange.min.max <- matrix(data=NA,ncol=2,nrow=ncol(PCmatrix))

        for(i in 1:ncol(PCmatrix)){
          PCrange.min.max[i,1] <- range(PCmatrix[,i])[1]
          PCrange.min.max[i,2] <- range(PCmatrix[,i])[2]
        }
        #new matrix for the difference between the max and min
        PCrange<-matrix(data=NA,ncol=1,nrow=ncol(PCmatrix))

          for (j in 1:ncol(PCmatrix)){
            PCrange[j,1] <- (PCrange.min.max[j,2]-(PCrange.min.max[j,1]))
          }
        return(PCrange)
  }

#------------------------------------------------------------
#Calculate the sum of ranges of each PC axis
  PCsumrange <- function(PCdata){
    ranges <- PCrange(PCdata)
    sumrange <- sum (ranges)
    return (sumrange)
  }

#------------------------------------------------------------
#Calculate the product of ranges of each PC axis
  PCprodrange <- function(PCdata){
    PCmatrix <- as.matrix(PCdata)
    ranges <- PCrange(PCmatrix)
      prodrange <- prod(ranges[,1])
        prodrange.scaled <- prodrange^(1/ncol(PCmatrix))
    return(prodrange.scaled)
  }   

################################################
#2) Based on interlandmark distances
################################################

#Claude (Morphometrics with R, 2008) function to calculate the interlandmark distances between pairs of coordinates 
  ild <- function(M1,M2){
         sqrt(sum((M1-M2)^2))
         }

#--------------------------------------------------------
#Calculate interlandmark distances between species' mean coordinates and reference shapes
  dist.to.ref <- function(coordinates){
    ref <- mshape(coordinates)
    ild.dist <- NULL
          for (i in 1:(dim(coordinates)[3])){
            ild.dist[i] <- ild(coordinates[,,i], ref)
          }
    return(ild.dist)
   }
   
#-----------------------------------------------------------------
#Calculate the morphological disparity; sum of squared distances/n-1, Zelditch 2012
  ZelditchMD <- function(distances){
                sum(((distances)^2)/(length(distances)-1))
                }      