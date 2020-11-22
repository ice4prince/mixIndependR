#'Calculate numbers of sharing alleles each pair at each locus
#'@details This function calculates the numbers of shared alleles between each pair of individuals for a dataset.Output a table and Usually followed by write.csv(as.data.frame(y),file = "~/*.csv") to export the results.
#'@importFrom stringr str_c
#'@usage AlleleShare_Table(x,replicate=TRUE)
#'@param x a dataset of alleles. Each row denotes each individual.One allele in one cell.In the (2r-1)th column, there is the same locus with the 2r-th column; noted: no column for ID, make row.names=1 when importing.
#'@param replicate a logical variable. if replicate is TRUE, the pairs are formed with replicates; if FALSE, the pairs are formed without replicate.
#'@return y a matrix of numbers of shared alleles. Each row denotes each pair; Each column denotes each locus.
#'@export
#'@examples
#'x <- data.frame(STR1=c(12,13,13,14,15,13,14,12,14,15),
#'                STR1_1=c(12,14,13,15,13,14,13,12,14,15),
#'                SNP1=c("A","T","A","A","T","A","A","T","T","A"),
#'                SNP1_1=c("A","T","T","T","A","T","A","A","T","T"))
#'AlleleShare_Table(x,replicate=TRUE)
#'

AlleleShare_Table <- function(x,replicate=TRUE)
{
  n <- nrow(x)
  m <- ncol(x)/2
  MakerTwo<-colnames(x)
  d <- mat.or.vec(2,n/2)
  if(replicate){
    d<-combn(n,2,simplify = TRUE)}else{
      if ((n %% 2) == 0){
        d[1,]<- sample(n,n/2,replace = FALSE)
        d[2,]<- sample(setdiff(c(1:n),d[1,]))
      }else{
        d[1,]<- sample(n,(n-1)/2,replace = FALSE)
        d[2,]<- sample(setdiff(c(1:n),d[1,]))[-1]
      }

    }

  p <- ncol(d) ###number of pairs####
  A.S. <- mat.or.vec(p,m)
  Name <- as.character(rep(NA, p))
  Namec <- rep.int(0,m)

  for (i in 1:p)
  {
    a <- as.character(x[d[1,i],])
    b <- as.character(x[d[2,i],])
    c <- rep.int(0,m)
    for (r in 1:m)
    {
      Namec[r]<-MakerTwo[(2*r-1)]
      if(identical(a[(2*r-1):(2*r)],b[(2*r-1):(2*r)])){
        c[r] <- 2
      }

      else{
        c[r]<-length(intersect(a[(2*r-1):(2*r)],b[(2*r-1):(2*r)]))
      }
    }
    A.S.[i,] <-c
    Name[i] <- str_c(d[1,i],d[2,i],sep = "vs")
  }

  colnames(A.S.) <- Namec
  row.names(A.S.)<- Name
  return(A.S.)
}
