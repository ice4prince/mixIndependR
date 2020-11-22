#'Generate a Bundle of Simulated distributions for No. of heterozygous loci with known heterozygosites
#'@details This function generates multinomial distribution for loci known the heterozygosity and build the simulated distribution for no. of heterozygous loci.
#'@importFrom stats chisq.test rbinom
#'@importFrom utils combn
#'@usage Simulate_DistK(H,m,t)
#'@param H a vector of average heterozygosity of each loci. Length of H is the number of loci.
#'@param m the sample size you want, usually similar to the real sample size.
#'@param t the number of samples you want to build
#'@return a matrix of frequencies of No. of Heterozygous Loci. Each row denotes each simulated sample; Each column denotes each No. of Heterozygous loci, from 0 to length of H.
#'@export
#'@examples
#'Simulate_DistK(runif(10),500,100)
#'

Simulate_DistK <- function(H,m,t){
  n <- length(H)
  OneSample <- mat.or.vec(n*t,m)
  Freq <- mat.or.vec(t,n+1) ####each row denotes each small sample, the column denotes the numbers of heterozygous loci: 0~n, the value in cell denotes how many individuals has such number of heterozygous loci.####
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
  for (i in 1:t){
    a <- ((i-1)*n)+1
    b <- i*n
    OneSample[(a:b),]<-replicate(m,rbinom(n,1,H))
    for (j in 0:n){
      Freq[i,j+1] <- countif(colSums(OneSample[(a:b),]),j)
    }
  }
  colnames(Freq) <- c(0:n)
  return(Freq)
}
