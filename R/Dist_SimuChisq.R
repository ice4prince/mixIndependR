#'Build a simulated distribution for Chi-Square
#'@details This function build the distribution of Chi square statistics for simulated samples
#'@usage Dist_SimuChisq(s,prob,b)
#'@param s a matrix of frequencies for each simulated sample. Each row for each sample.
#'@param prob a vector of expected probability for each simulated sample.
#'@param b the times of bootstrapping.
#'@return a vector of Chi-square statistics, length is the times of sampling.
#'@export
#'@examples
#'require(mixIndependR)
#'h<-runif(10)
#'s<-Simulate_DistK(h,500,100)
#'Exp <- DistHetero(h)
#'Dist_SimuChisq(s,Exp$Density,10)

Dist_SimuChisq <- function(s,prob,b){
  idx <-which(prob==0)
  if (length(idx)==0){
    prob <- prob
    s <- s
  }else{
    prob <- prob[-idx]
    s <- s[,-idx]
  }
  t <- nrow(s) ####times of sampling.######
  MyChisqValues <- replicate(t,0)
  for (k in 1:t){
    MyChisqValues[k] <- chisq.test(s[k,],p=prob,simulate.p.value = T,B=b)$statistic
  }
  return(MyChisqValues)
}
