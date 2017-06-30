#' Compute the correlation matrix based on a PowExp structure
#' 
#' The function computes the m by n correlation between points in the hold set and points in the training set. 
#' 
#' @param dis.tem a m*n by d matrix contains the m*n combinations of (x_hk - x_jk) for h <= m, j <=n, and k < d. 
#' Must be passed from cov.r1.dis.power.full().   
#' @param lambda sensitivity parameters at lambda scale.
#' @param n size of training set. 
#' @param m size of test set. 
#' @param p smoothness parameters from 1 to 2.
#' @return r.matrix returns the correlation matrix between all points in the test set and all points in the training set.
#'  

cov.r1.power.full<-function(dis.tem,lambda,p, n, m){
  rho=exp(lambda)/(1+exp(lambda))
  phi=-4*log(rho); 
  dis<-NULL
  for(i in 1:dim(dis.tem)[2]){
    dis<-cbind(dis, (abs((dis.tem[,i])))^(p[i]))  
  }  
  r<-exp(-(dis%*%phi))
  r.matrix<-matrix(r,nrow=m,ncol=n,byrow=T)
  return(r.matrix)
}
