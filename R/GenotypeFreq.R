#'Calculate Genotype Frequency###
#'@details This function calculates the observed or expected genotype frequency from dataset and allele frequency.#####
#'@usage GenotypeFreq(x,sep,expect=TRUE)
#'@param x a dataframe of genotype data with rownames of sample ID and column names of markers.
#'@param sep allele separator in the imported genotype data. Note: when using the special character like "|", remember to protect it as "\\|"(default).
#'@param expect a logic variable. If expect is true, the function will calculate the expected genotype probabilities. If false, calculate the observed genotype frequencies.
#'@references Chakraborty, R., Srinivasan, M. R., & Daiger, S. P. (1993, ISSN:0002-9297).
#'@return a matrix of genotype frequencies. Each row denotes each genotype; each column denotes each loci. The order of markers follows x; the genotypes are ordered by: from 1:l-th column, the genotypes are homozygous in order as : p1p1, p2p2,p3p3,...,plpl;from ll-th to u-th column, the genotypes are heterozygous in order as:choose(l,2) like: p1p2,p1p3,...,p1pl,p2p3,p2p4,...p2pl,...p(l-1)pl
#'@export
#'@examples
#'require(mixIndependR)
#'x <- data.frame(STR1=c("12|12","13|14","13|13","14|15","15|13","13|14","14|13","12|12","14|14","15|15"),
#'                 SNP1=c("A|A","T|T","A|T","A|T","T|A","A|T","A|A","T|A","T|T","A|T"))
#'GenotypeFreq(x,"\\|",expect=TRUE)
#'
#'
#'
GenotypeFreq <- function(x,sep="\\|",expect = TRUE){
  p <-AlleleFreq(x,sep)
  Gt_a<-as.data.frame(cbind(rbind(rownames(p),rownames(p)),combn(rownames(p),2)))
  Gt<-as.vector(sapply(Gt_a,function(x){paste0(x[1],"|",x[2])}))
  if(expect){
    ho<-p*p
    he<-sapply(data.frame(p),function(x){combn(x,2)[1,]*combn(x,2)[2,]})
    output <-rbind(ho,he)
  }else{
    output <-sapply(x,function(x){sapply(Gt,counta,z=x)})
  }
  rownames(output) <- Gt
  return(output)
}
