#'Quick pvalue of total number of heterozygous loci
#'@details This function is a summary of pipeline for number of heterozygous loci (K), and generates the p-value of K for the target dataset.
#'@usage mixIndependK(x,t,B)
#'@importFrom stats ecdf
#'@param x a dataset of alleles. Each row denotes each individual.One allele in one cell.In the (2r-1)th column, there is the same locus with the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@param t times of simulation in "Simulate_DistK" and "Simulate_DistX".
#'@param B times of bootstrapping in Chi Squares Test.
#'@return pvalue (1-cumulative probabilities) for the number of heterozygous loci(K)
#'@export
#'@examples
#'x0 <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                 STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                 SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                 SNP1_1=c("A","T","T","T","A","T","A","A","T","T"),
#'                 STR2=c(10,12,11,9,10,12,11,12,12,10),
#'                 STR2_1=c(10,9,11,11,10,12,10,10,12,9),
#'                 SNP2=c("C","C","G","G","G","G","C","G","G","C"),
#'                 SNP2_1=c("C","C","G","G","C","G","C","C","G","G"))
#'mixIndependK(x0,10,10)

mixIndependK<-function(x,t,B){
  s <- nrow(x)
  p <- AlleleFreq(x)
  h <- Heterozygous(x)
  H <- RxpHetero(h,p,HWE = F)
  Obs_DistHetero<-FreqHetero(h)
  Exp_DistHetero<-DistHetero(H)
  prob<-Exp_DistHetero$Density
  obs<-Obs_DistHetero$Frequency
  idx <-which(prob==0)
  if (length(idx)==0){
    prob <- prob
    obs <- obs
  }else{
    prob <- prob[-idx]
    obs <- obs[-idx]
  }
  x20 <-chisq.test(obs,p=prob,simulate.p.value = T,B=B)
  s<-Simulate_DistK(H,s,t)
  x2<-Dist_SimuChisq(s,Exp_DistHetero$Density,B)
  P <- ecdf(x2)
  return(1-P(x20$statistic))
}
