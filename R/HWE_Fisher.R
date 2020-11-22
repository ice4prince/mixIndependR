#'Test the Hardy Weinberg Equilibrium with Fisher's exact test###
#'@details This function check the Hardy Weinberg Equilibrium with Fisher's exact Test.#####
#'@usage HWE.Fisher(p,H,y)
#'@param p a matrix of allele frequency;each row denotes allele; each column denotes each loci;
#'@param H a vector of number of Heterozygotes on each loci; length is number of loci.
#'@param y a matrix of observed genotype Densities(Not count). Each row denotes each genotype; each column denotes each loci. The order of markers follows x; the genotypes are ordered by: from 1:l-th column, the genotypes are homozygous in order as : p1p1, p2p2,p3p3,...,plpl;from ll-th to u-th column, the genotypes are heterozygous in order as:choose(l,2) like: p1p2,p1p3,...,p1pl,p2p3,p2p4,...p2pl,...p(l-1)pl
#'@return a vector of p-values of Fisher's test; ordered by the order of loci in p or x
#'@references Weir, B. S. (1996, ISBN:9780878939022)
#'@export
#'@examples
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'p <-AlleleFreq(x)
#'G <- GenotypeFreq(x,p,expect = FALSE)
#'h <- Heterozygous(x)
#'H <- RxpHetero(h,p,HWE = FALSE)
#'HWE.Fisher(p,H,G/colSums(G))
#'
HWE.Fisher <- function(p,H,y){
  n<- colSums(y)[1] #####number of samples#####
  m <- ncol(y) ####number of loci#####
  p0 <- p*2*n
  Pr <- rep.int(0,m)
  for (i in 1:m){
    Pr[i]<-exp(lfactorial(n)+H[i]*log(2)+sum(lfactorial(p0[,i]))-lfactorial(2*n)-sum(lfactorial(y[,i])))
  }
  names(Pr) <- colnames(y)
  return(Pr)

}
