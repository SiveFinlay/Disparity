
#09/05/2014 

#Updated original functions code (29/01/2014) 
  #Made the PCvariance functions easier to read
  #Added the ild.distance and ZelditchMD functions to calculate disparity based on interlandmark distances (Zelditch 2012)
  
#Original code for disparity functions in the DisparityPractice_16_01_14_using_SkLat_08_11_13 script

#PCvariance: calculates the variance of PC axes
#PCsumvar: sum of variance of PC axes
#PCprodvar: product of variance of PC axes

#PCrange: range of PC axes
#PCsumrange: sum of the range of PC axes
#PCprodrange: product of the range of PC axes

#ild: interlandmark distance (from Claude 2008)
#dist.to.ref: ild between species' mean coorinates and the overall mean (reference) shape
#ZelditchMD: disparity as sum of squared distances/n-1


#################################################
#NB: I checked all of these functions by comparing their output to my manual calculation of the
#sum and product of range and variance for tenPC axes in the other script file (DisparityPractice_using_SkLat_08_11_13)
#They all gave the same answers so these functions are working in the way that I want them to
#############################################################
#Variance

#the first option is to have a function which just calculates the variances of each of the PC axes
#where x1 is a data matrix (number of rows is the number of individuals and number of columns is the number of PC axes used)
PCvariance<-function(PCmatrix){
  cols<-ncol(PCmatrix)
  if (class(PCmatrix)=="matrix"){
    PCvar<-matrix(data=NA,ncol=cols,nrow=1)
      for (i in 1:cols){
        PCvar[,i]<-var(PCmatrix[,i])
      }
  return(PCvar)
  }
  else (stop("This function requires data of class type matrix"))
  }

#then just get the sum and product of the variance in a separate step

#second option is a function specifically to calculate the sum of variance; returns one output
#(the cat options return two outputs but then I can't select them from the result)
PCsumvar<-function(PCmatrix){
  cols<-ncol(PCmatrix)
  if (class(PCmatrix)=="matrix"){
    PCvar<-matrix(data=NA,ncol=cols,nrow=1)
      for (j in 1:cols){
        PCvar[,j]<-var(PCmatrix[,j])
      }
    sumvar<-sum(PCvar)
    return(sumvar)
  }
  else (stop("This function requires data of class type matrix"))
  }

#and the final variance function is one to calculate the product of the variance directly
PCprodvar<-function(PCmatrix){
  cols<-ncol(PCmatrix)
  if (class(PCmatrix)=="matrix"){
    PCvar<-matrix(data=NA,ncol=cols,nrow=1)
      for (j in 1:cols){
        PCvar[,j]<-var(PCmatrix[,j])
      }
    prodvar<-prod(PCvar)
    #the Wills and Brusatte papers scale their product values
    #divide the product by the root of the number of axes used to reduce the dimensionality of the answer
      prodvar.scaled<-prodvar^(1/cols)
    #return the scaled product since that's the useful one
    return(prodvar.scaled)
    }
  else (stop("This function requires data of class type matrix"))
  }

#########################################################
#Functions for calculating the range, sum of ranges and product of ranges
#same as above, x1 is a matrix where the number of rows is the number of individuals and the
#number of columns is the number of PC axes used)

#Function to just calculate the range
PCrange<-function(PCmatrix){
  cols<-ncol(PCmatrix)
  if (class(PCmatrix)=="matrix"){
    #empty matrix for the minimum and maximum values of each PC axis
    PCrange.min.max<-matrix(data=NA,ncol=2,nrow=cols)
      for(i in 1:cols){
        #put the minimum value into the first column
        PCrange.min.max[i,1]<-range(PCmatrix[,i])[1]
        #and the maximum value into the second column
        PCrange.min.max[i,2]<-range(PCmatrix[,i])[2]
      }
    #make a new matrix for the single range values (difference between the max and min)
    PCrange<-matrix(data=NA,ncol=1,nrow=cols)
      #fill the matrix with the range (difference between min and max)
      for (j in 1:cols){
        PCrange[j,1]<-(PCrange.min.max[j,2]-(PCrange.min.max[j,1]))
        #the minimum values are often negative so I need to put them into brackets for the
        #subtraction from the maximum value
        #otherwise it doesn't calculate the full range (difference between 3-2 and 3-(-2))
      }
    return(PCrange)
  }
  else (stop("This function requires data of class type matrix"))
}

#I could just use the output of this function to get the sum and products of the ranges separately
#Alternatively, here are two functions to get the sum and products of the ranges directly

PCsumrange<-function(PCmatrix){
  cols<-ncol(PCmatrix)
  if (class(PCmatrix)=="matrix"){
  #empty matrix for the minimum and maximum values of each PC axis
    PCrange.min.max<-matrix(data=NA,ncol=2,nrow=cols)
        for(i in 1:cols){
          #put the minimum value into the first column
          PCrange.min.max[i,1]<-range(PCmatrix[,i])[1]
          #and the maximum value into the second column
          PCrange.min.max[i,2]<-range(PCmatrix[,i])[2]
        }
    #make a new matrix for the single range values (difference between the max and min)
    PCrange<-matrix(data=NA,ncol=1,nrow=cols)
    #fill the matrix with the range (difference between min and max)
      for (j in 1:cols){
        PCrange[j,1]<-(PCrange.min.max[j,2]-(PCrange.min.max[j,1]))
      }
    sumrange<-sum(PCrange[,1])
    return(sumrange)
  }
  else (stop("This function requires data of class type matrix"))
}


#and the same idea for the product
PCprodrange<-function(PCmatrix){
  cols<-ncol(PCmatrix)
    if (class(PCmatrix)=="matrix"){
    #empty matrix for the minimum and maximum values of each PC axis
      PCrange.min.max<-matrix(data=NA,ncol=2,nrow=cols)
      for(i in 1:cols){
        #put the minimum value into the first column
         PCrange.min.max[i,1]<-range(PCmatrix[,i])[1]
        #and the maximum value into the second column
        PCrange.min.max[i,2]<-range(PCmatrix[,i])[2]
      }
    #make a new matrix for the single range values (difference between the max and min)
    PCrange<-matrix(data=NA,ncol=1,nrow=cols)
    #fill the matrix with the range (difference between min and max)
    for (j in 1:cols){
      PCrange[j,1]<-(PCrange.min.max[j,2]-(PCrange.min.max[j,1]))
    }
    prodrange<-prod(PCrange[,1])
    #scaled product (divide by the root of the number of PC axes)
    prodrange.scaled<-prodrange^(1/cols)
    #return the scaled product since that's the useful one
    return(prodrange.scaled)
  }
  else (stop("This function requires data of class type matrix"))
  }


#-----------------------------------------------
#Claude (Morphometrics with R, 2008) function to calculate the interlandmark distances between pairs of coordinates 
  ild<-function(M1,M2){
        sqrt(sum((M1-M2)^2))
        }


#Function to calculate interlandmark distances between species' mean coordinates and reference shapes

dist.to.ref<-function(coordinates, fam, species){
  ref<-mshape(coordinates)
   ild.dist<-as.data.frame(matrix(NA, ncol=3, nrow=(length(species))))
      colnames(ild.dist)<-c("Family","Binomial","Ild")
      ild.dist[,1]<-fam
      ild.dist[,2]<-species
        for (i in 1:length(binom)){
              ild.dist[i,3]<-ild(coordinates[,,i], ref)
       }
    return(ild.dist)
   }
   

#Function to calculate the morphological disparity; sum of squared distances/n-1
ZelditchMD<-function(distances){
                sum(((distances)^2)/(length(distances)-1))
                }      