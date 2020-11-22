#'Calculate the Real Probability of 0,1 and 2 Shared Alleles###
#'@details This function Calculates the density of 0,1 and 2 Shared Alleles for a set of loci. Usually followed by write.csv(as.data.frame(y),file = "~/*.csv") to export the result of a n x3 matrix.
#'@usage RealProAlleleShare(AS)
#'@param AS a matrix/double of no. of Shared alleles, made up with 0,1 and 2; Outcome of "AlleleShare_Table". Each column denotes each locus. Each row denotes each individual.
#'@return a matrix/double of real density of 0,1 and 2 shared alleles for each locus. Each row denotes each locus. The first column denotes the probability of 0 shared alleles, the second denotes 1 shared allele, the third denotes 2 shared alleles.
#'@export
#'@examples
#'AS<-matrix(sample(c(0:2),20,replace=TRUE,prob=c(0.3,0.3,0.4)),nrow=5)
#'RealProAlleleShare(AS)
RealProAlleleShare <-function(AS){
  l <-nrow(AS)
  n <-ncol(AS)
  Pro_real <- mat.or.vec(n,3)
  for (i in 1:n){
    Pro_real[i,1]<-counta(AS[,i],0)
    Pro_real[i,2]<-counta(AS[,i],1)
    Pro_real[i,3]<-counta(AS[,i],2)
  }
  row.names(Pro_real) <- colnames(AS)
  colnames(Pro_real) <- c("P0","P1","P2")
  return(Pro_real/l)
}
