#'Test the Hardy Weinberg Equilibrium with Chi-square test####
#'@details This function check the Hardy Weinberg Equilibrium from observed and expected distribution with Chi-square test#####
#'@importFrom stats chisq.test
#'@usage HWE.Chisq(x,x0,rescale.p=FALSE,simulate.p.value=FALSE,B)
#'@param x a matrix of observed genotype frequencies. Each row denotes each genotype; each column denotes each loci. The order of markers follows x; the genotypes are ordered by: from 1:l-th column, the genotypes are homozygous in order as : p1p1, p2p2,p3p3,...,plpl;from ll-th to u-th column, the genotypes are heterozygous in order as:choose(l,2) like: p1p2,p1p3,...,p1pl,p2p3,p2p4,...p2pl,...p(l-1)pl
#'@param x0 a matrix of expected Probabilities;each row denotes each genotype; each column denotes each loci. The order of markers follows x; the genotypes are ordered by: from 1:l-th column, the genotypes are homozygous in order as : p1p1, p2p2,p3p3,...,plpl;from ll-th to u-th column, the genotypes are heterozygous in order as:choose(l,2) like: p1p2,p1p3,...,p1pl,p2p3,p2p4,...p2pl,...p(l-1)pl
#'@param rescale.p a logical scalar; if TRUE then p is rescaled (if necessary) to sum to 1. If rescale.p is FALSE, and p does not sum to 1, an error is given.
#'@param simulate.p.value a logical indicating whether to compute p-values by Monte Carlo simulation.
#'@param B an integer specifying the number of replicates used in the Monte Carlo test.
#'@return y a list of result of chi-square test, $chi $pvalue; chi and pvalue are vectors of chi square statistics/ p values. Orders follows x.
#'@export
#'@examples
#'require(mixIndependR)
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'p <- AlleleFreq(x)
#'g <- GenotypeFreq(x,p,expect=FALSE)
#'g0 <- GenotypeFreq(x,p,expect=TRUE)
#'HWE.Chisq(g,g0,rescale.p=FALSE,simulate.p.value=TRUE,2000)

HWE.Chisq <- function(x,x0,rescale.p=FALSE,simulate.p.value=FALSE,B){
  n <- ncol(x)
  chi <- rep.int(0,n)
  pvalue <- rep.int(0,n)
  if(rescale.p){
    if (simulate.p.value){
      for (i in 1:n){
        ss <-x[,i]
        so <-x0[,i]
        cr <- 5/sum(ss)
        if(length(ss[which(so>cr)])>2){
          chi[i]<-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = T,simulate.p.value = T, B=B)$statistic
          pvalue[i] <-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = T,simulate.p.value = T, B=B)$p.value
        }else{
          chi[i]<-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = T,simulate.p.value = T, B=B)$statistic
          pvalue[i] <-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = T,simulate.p.value = T, B=B)$p.value
        }
              }
    }else{
      for (i in 1:n){
        ss <-x[,i]
        so <-x0[,i]
        cr <- 5/sum(ss)
        if(length(ss[which(so>cr)])>2){
          chi[i]<-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = T,simulate.p.value = F, B=B)$statistic
          pvalue[i] <-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = T,simulate.p.value = F, B=B)$p.value
        }else{
          chi[i]<-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = T,simulate.p.value = F, B=B)$statistic
          pvalue[i] <-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = T,simulate.p.value = F, B=B)$p.value
        }
        }
  }
    }else{
    if(simulate.p.value){
    for (i in 1:n){
      ss <-x[,i]
      so <-x0[,i]
      cr <- 5/sum(ss)
      if(length(ss[which(so>cr)])>2){
        chi[i]<-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = F,simulate.p.value = T, B=B)$statistic
        pvalue[i] <-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = F,simulate.p.value = T, B=B)$p.value
      }else{
        chi[i]<-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = F,simulate.p.value = T, B=B)$statistic
        pvalue[i] <-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = F,simulate.p.value = T, B=B)$p.value
      }
      }
      }else{
    for (i in 1:n){
      ss <-x[,i]
      so <-x0[,i]
      cr <- 5/sum(ss)
      if(length(ss[which(so>cr)])>2){
        chi[i]<-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = F,simulate.p.value = F, B=B)$statistic
        pvalue[i] <-chisq.test(ss[which(so>cr)],p=so[which(so>cr)]/sum(so[which(so>cr)]),rescale.p = F,simulate.p.value = F, B=B)$p.value
      }else{
        chi[i]<-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = F,simulate.p.value = F, B=B)$statistic
        pvalue[i] <-chisq.test(ss[which(ss>0)],p=so[which(ss>0)]/sum(so[which(ss>0)]),rescale.p = F,simulate.p.value = F, B=B)$p.value
      }
      }
      }
  }
  names(chi) <- colnames(x)
  names(pvalue) <- colnames(x)
  y <- list("chi"=chi,"pvalue"=pvalue)
  return(y)
  print(y)
  }
