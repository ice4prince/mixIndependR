#'Calculate Genotype Frequency###
#'@details This function calculates the observed or expected genotype frequency from dataset and allele frequency.#####
#'@usage GenotypeFreq(x,p,expect=TRUE)
#'@param x a dataset of alleles. Each row denotes each sample. One allele in one cell.In the (2r-1)th column, there is the other allele on the same locus from that in the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@param p a matrix of allele frequencies. Each row denotes each allele; each column denotes each marker. The order of markers follows x.
#'@param expect a logic variable. If expect is true, the function will calculate the expected genotype probabilities. If false, calculate the observed genotype frequencies.
#'@references Chakraborty, R., Srinivasan, M. R., & Daiger, S. P. (1993, ISSN:0002-9297).
#'@return y a matrix of genotype frequencies. Each row denotes each genotype; each column denotes each loci. The order of markers follows x; the genotypes are ordered by: from 1:l-th column, the genotypes are homozygous in order as : p1p1, p2p2,p3p3,...,plpl;from ll-th to u-th column, the genotypes are heterozygous in order as:choose(l,2) like: p1p2,p1p3,...,p1pl,p2p3,p2p4,...p2pl,...p(l-1)pl
#'@export
#'@examples
#'require(mixIndependR)
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'p <- AlleleFreq(x)
#'GenotypeFreq(x,p,expect=TRUE)
#'
#'
#'
GenotypeFreq <- function(x,p,expect = TRUE){
  x <- as.matrix(x)
  n <- ncol(p) ###number of loci;
  l <- nrow(p)  ###number of alleles;
  m <- nrow(x)  ####number of individuals/Samples;
  ll <- l+1
  c <- choose(l,2)+l   ####number of heterozygous and homozygous genotypes###
  cb <-combn(c(1:l),2)
  d<-row.names(p)
  het<-combn(d,2)
  homo <-matrix(c(combn(d,1),combn(d,1)),nrow = 2,byrow = T)
  Gn<-cbind(homo,het)
  y <- mat.or.vec(c,n)
  if(expect){
    for (i in 1:n){
      for (h in 1:l){
        y[h,i] <- p[h,i]^2
      }
      for (h in ll:c){
        y[h,i] <- 2*p[cb[1,(h-l)],i]*p[cb[2,(h-l)],i]
      }
    }
  }
  if(!expect){
    for (k in 1:n){
      a <- 2*k -1
      b <- 2*k
      for (j in 1:c){
        y[j,k] <- 0
        for (r in 1:m){
          if (setequal(x[r,a:b],Gn[,j])){
            y[j,k] <- y[j,k]+1
          }else{
            y[j,k] <- y[j,k]
          }
        }
      }
    }
  }
  colnames(y) <- colnames(p)
  s <- ncol(Gn)
  G.Nam <- rep.int(s,0)
  for (g in 1:s){
    G.Nam[g]<- paste(Gn[1,g],Gn[2,g],sep = " ")
  }
  rownames(y) <- G.Nam
  return(y)
}
