#'Build Observed Distribution of No. of Heterozygous loci
#'@details This function build the observed distributions from observed heterozygosity table, made up of 0,1.
#'@usage FreqHetero(h)
#'@param h a dataframe of heterozygosity, made up with 0 and 1, outcome of function "Heterozygous" Rows for individuals, and columns for markers.
#'@return a dataframe of frequencies of each number of heterozygous loci(from 0 to No. of loci)
#'@export
#'@examples
#'h<-matrix(rbinom(20,1,0.5),nrow=5)
#'FreqHetero(h)

FreqHetero<- function(h){
  n<-ncol(h)
  K<-rowSums(h)
  Frequency <- rep.int(0,(n+1))
  for (i in 1:(n+1)){
    Frequency[i] <- counta(K,(i-1))
  }
  K_io <- c(0:n)
  df <- data.frame(K_io,Frequency)
  return(df)
}
