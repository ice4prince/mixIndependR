#'Build a simulated distribution for No. of Shared Alleles
#'@details This function generates multinomial distribution for loci known the Allele Frequency and Expected Probability of Shared 2,1 or 0 alleles
#'@usage Simulate_DistX(e,m,t)
#'@param e a matrix of Probability of Sharing 2,1 or 0 alleles at each loci. Each row denotes each locus. Three columns denote sharing 0,1 or 2 alleles.
#'@param m the sample size you want, usually similar to the real sample size.
#'@param t the number of samples you want to build/ the times to generate a sample
#'@return y a matrix of frequencies of No. of shared alleles. Each row denotes each simulated sample; Each column denotes each No. of shared alleles, from 0 to 2e length of e.
#'@export
#'@examples
#'e0<-data.frame("P0"=runif(5,min = 0,max = 0.5),"P1"=runif(5,0,0.5))
#'e<-data.frame(e0,"P2"=1-rowSums(e0))
#'Simulate_DistX(e,500,10)
#'

Simulate_DistX <- function(e,m,t){
  n <- nrow(e)
  OneSample <- mat.or.vec(n*t,m)
  Freq <- mat.or.vec(t,2*n+1) ####each row denotes each small sample, the column denotes the numbers of shared alleles: 0~2n, the value in cell denotes how many individuals has such number of shared alleles.####
  ###build the function for count##
  countif <- function(z,y){
    conts <- 0
    s <- length(z)
    for (l in 1:s){
      if (z[l]==y){
        conts <-conts+1
      }
    }
    return(conts)
  }
  #####pick up m individuals t times, to make up the t small samples. ###
  for (j in 1:t){
    for (i in 1:n){
      OneSample[((j-1)*n+i),] <- replicate(m,sample(c(0,1,2),1,replace=T,prob = e[i,]))
    }
    a <- ((j-1)*n)+1
    b <- j*n
    N <- 2*n
    for (k in 0:N){
      Freq[j,k+1] <- countif(colSums(OneSample[(a:b),]),k)
    }
  }
  colnames(Freq) <- c(0:(2*n))
  return(Freq)
}
