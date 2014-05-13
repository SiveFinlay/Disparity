#29/01/14
#Function to calculate a p value from a distribution
#Copied this code from the DisparityPractice_16_01_14_using_SkLat_08_11_13 script

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
#--------------------------------------------------------
#test the function
#make a random normal distribution to test out a function
#  rand.dist<-(rnorm(1:100, mean=5, sd=1))
#  range(rand.dist)
#test with an observed value at the lower end of the range
#  pvalue.dist(distribution=rand.dist,obs.val=3)
#test with an observed value at the higher end of the range
#  pvalue.dist(distribution=rand.dist,obs.val=6.5)
#each time the function selects the correct option of lowerp vs higherp as the right output value