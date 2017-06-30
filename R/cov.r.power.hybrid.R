#' Compute the correlation matrix based on a PowExp structure
#' 
#' The function computes the m by n correlation between points in the hold set and points in the training set. 
#' 
#' @param dis a m*n by d matrix contains the m*n combinations of (x_hk - x_jk)^p_k for h <= m, j <=n, and k < d. 
#' Must be passed from cov.r1.dis.power.hybrid().   
#' @param lambda sensitivity parameters at lambda scale.
#' @param n size of training set. 
#' @param m size of test set. 
#' @return r.matrix returns the correlation matrix between all points in the test set and all points in the training set.
#'  

cov.r1.power.hybrid <- function(dis, lambda, n, m){
  rho=exp(lambda)/(1+exp(lambda))
  phi=-4*log(rho); 
  dis<-dis
  r<-exp(-(dis%*%phi))
  r.matrix<-matrix(r,nrow=m,ncol=n,byrow=T)
  return(r.matrix)
}