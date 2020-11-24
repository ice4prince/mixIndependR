#'Split Genotype Table to Duo-Allele Table
#'@description Split each column to two columns for a table of genotypes
#'@usage split_Genotype(df,sep="\\|",dif="*",byrow=TRUE)
#'@import tidyr
#'@param df a dataframe of genotype data with rownames of sample ID and column names of markers.
#'@param sep allele separator in the imported genotype data. Note: when using the special character like "|", remember to protect it as "\\|"(default).
#'@param dif a symbol differentiate the one marker on each allele.
#'@param byrow a logical variable. If byrow is TRUE, the output is arranged with double rows but the same columns, and the table of the second allele is followed after the first allele table by rows with double individual IDs in the same order.\cr If byrow is false,the output is arranged by double columns and the same rows; the column names are in the order of alphabet by pairs.
#'@details The function convert a genotype data to allele data with double columns or with double rows; the rownames are sample ID in the same order but twice if the rows are doubled, and the column names are in the same order or in the order of alphabet by pairs if columns are doubled. \cr The paremater "sep" is the symbol of allele separator in the imported genotype data. \cr The parameter "dif" is the difference between the second and the first apprearance for the same marker. For example, if "dif = _1", the column names of output will be "marker1" "marker1 _1","marker2","marker2 _1", if the original list of column names is "marker1","marker2".
#'@return a dataframe with doubled columns of import data and alleles in different columns
#'@examples
#'df <- data.frame(SNP1=c("A|A","T|T","A|T","A|T"),
#'                 STR1=c("12|12","13|14","13|13","14|15"))
#'split_Genotype(df,sep="\\|",dif="*",byrow=TRUE)

split_Genotype <- function(df,sep="\\|",dif="*",byrow=TRUE){
  g0 <- as.matrix(df)
  n <- nrow(g0)
  m <- ncol(g0)
  a1<-matrix(sapply(strsplit(g0,sep),"[",1),nrow = n,ncol=m,byrow = F)
  a2<-matrix(sapply(strsplit(g0,sep),"[",2),nrow = n,ncol=m,byrow = F)
  colnames(a1) <- colnames(df)
  rownames(a1) <- rownames(df)
  colnames(a2) <- paste(colnames(df),dif)
  rownames(a2) <- paste(rownames(df),dif)
  if (byrow){
    a <- rbind(a1,a2)
    output <-as.data.frame(a)
  }else{
    a <- cbind(a1,a2)
    output <-as.data.frame(a[,sort(colnames(a))])
  }
  return(output)
}
