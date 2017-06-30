#' Compute the predicted mean and variance based on GP with a PowExp structure
#' 
#' The function returns the predicted mean and variance.  
#' 
#' @param y True outputs of x. Size is n. 
#' @param lambda sensitivity parameters at lambda scale.
#' @param d a n by n matrix with all elements 0. Must be passed from dist.R.power.hybrid(). 
#' @param odut the indices by which to place the n*(n-1)/2 correlation 
#' in the correlation matrix filled by rows. Must be passed from dist.R.power.hybrid()  
#' @param dup a n*(n-1)/2 by d matrix contains the n*(n-1)/2 combinations 
#' of (x_ik - x_jk)^p_k for any i < j and k < d. Must be passed from dist.R.power.hybrid() 
#' @param dis a m*n by d matrix contains the m*n combinations of (x_hk - x_jk)^p_k for any for h <= m, j <=n, and k < d.   
#' Must be passed from cov.r1.dis.power.hybrid(). 
#' @return res returns the predicted mean. 
#' @return v.term2 return the predicted variance.

pred1.power.hybrid <- function(d, odut, dup, dis, y, lambda){  
  rho=exp(lambda)/(1+exp(lambda))
  phi=-4*log(rho); 
  n<-dim(d)[1]
  m<-dim(dis)[1]/n
  one<-as.matrix(rep(1,n))
  y.vector<-as.matrix(y)
  res<-matrix(nrow=1,ncol=m)
  v.term2<-matrix(nrow=1,ncol=m)
  gpost<-log.post.power.hybrid(d, odut,dup, y, lambda)
  mu.hat<-as.numeric(gpost$mu)
  sigma.hat<-as.numeric(gpost$sigma)
  hat.R<-gpost$R
  rm(gpost)
  U<-chol(hat.R)
  r.matrix<-cov.r1.power.hybrid(dis, lambda, n, m)
  r.t<-t(r.matrix)
  cRr<-backsolve(U,r.t,transpose=T)
  cRyb<-backsolve(U,(y.vector-one*mu.hat), transpose=T)
  res<-t(matrix(rep(mu.hat,m),ncol=1)+crossprod(cRr,cRyb))
  v.term2<-t(s.square1(sigma.hat,r.matrix,hat.R))
  dd<-list(res=res, v.term2=v.term2)
  return(dd)
}
