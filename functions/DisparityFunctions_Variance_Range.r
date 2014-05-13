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
  PCvariance <- function(PCmatrix){
    cols<-ncol(PCmatrix)
      PCvar <- matrix(data=NA,ncol=cols,nrow=1)

        for (i in 1:cols){
          PCvar[,i] <- var(PCmatrix[,i])
        }
  return(PCvar)
  }

#------------------------------------------------------------
#Caculate the sum of variance
  PCsumvar <- function(PCmatrix){
    cols <- ncol(PCmatrix)
      PCvar <- matrix(data=NA,ncol=cols,nrow=1)

        for (j in 1:cols){
          PCvar[,j] <-var (PCmatrix[,j])
        }
      sumvar<-sum(PCvar)
      return(sumvar)
  }

#------------------------------------------------------------
#Calculate the product of variance
  PCprodvar <- function(PCmatrix){
    cols <- ncol(PCmatrix)
      PCvar <- matrix(data=NA,ncol=cols,nrow=1)

        for (j in 1:cols){
          PCvar[,j] <- var(PCmatrix[,j])
        }
      prodvar<-prod(PCvar)
      #divide the product by the root of the number of axes used to reduce the dimensionality of the answer (cf Wills + Brusatte)
        prodvar.scaled <- prodvar^(1/cols)
        return(prodvar.scaled)
  }

#------------------------------------------------------------
#Calculate the range of each PC axis
  PCrange <- function(PCmatrix){
    cols <- ncol(PCmatrix)
      PCrange.min.max <- matrix(data=NA,ncol=2,nrow=cols)

        for(i in 1:cols){
          PCrange.min.max[i,1] <- range(PCmatrix[,i])[1]
          PCrange.min.max[i,2] <- range(PCmatrix[,i])[2]
        }
        #new matrix for the difference between the max and min
        PCrange<-matrix(data=NA,ncol=1,nrow=cols)

          for (j in 1:cols){
            PCrange[j,1] <- (PCrange.min.max[j,2]-(PCrange.min.max[j,1]))
          }
        return(PCrange)
  }

#------------------------------------------------------------
#Calculate the sum of ranges of each PC axis
  PCsumrange <- function(PCmatrix){
    ranges <- PCrange(PCmatrix)
    sumrange <- sum (ranges[,1])
    return (sumrange)
  }

#------------------------------------------------------------
#Calculate the product of ranges of each PC axis
  PCprodrange <- function(PCmatrix){
    cols <- ncol(PCmatrix)
    ranges <- PCrange(PCmatrix)
      prodrange <- prod(ranges[,1])
        prodrange.scaled <- prodrange^(1/cols)
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
  dist.to.ref <- function(coordinates, fam, species){
    ref <- mshape(coordinates)
    ild.dist <- as.data.frame(matrix(NA, ncol=3, nrow=(length(species))))
      colnames(ild.dist) <- c("Family","Binomial","Ild")
        ild.dist[,1] <- fam
        ild.dist[,2] <- species

          for (i in 1:length(binom)){
            ild.dist[i,3] <- ild(coordinates[,,i], ref)
          }
    return(ild.dist)
   }
   
#-----------------------------------------------------------------
#Calculate the morphological disparity; sum of squared distances/n-1, Zelditch 2012
  ZelditchMD <- function(distances){
                sum(((distances)^2)/(length(distances)-1))
                }      