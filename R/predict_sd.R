#' Predicted variance of a GP
#' 
#' The function computes the predicted variances. 
#'  
#' @param r The m by n correlation matrix of all points in the training set and all points in the test set.
#' Must be passed from cov.r1.power.hybrid().    
#' @param sigma The estimated variance of the GP. Must be passed from log.post.power.hybrid().   
#' @param R The correlation matrix of the training set. Must be passed from log.post.power.hybrid(). 
#' @return rett returns the predicted variances of all points in the test set. 
#'

s.square1 <- function(sigma,r,R){
  sigma<-sigma;R<-R;r=r
  n=dim(R)[1]
  #m<-dim(x)[1]; 
  one<-as.matrix(rep(1,n),ncol=1)
  U<-chol(R)
  r.t<-t(r)
  Rinf=backsolve(U,r.t,transpose=T)
  #fRinf=crossprod(Rinf)
  #tem1<-1-as.numeric(diag(fRinf))
  #rm(fRinf)
  ###### changed to optimize the speed
  #tem1 = 1 - as.numeric(t(Rinf^2)%*%matrix(1, nrow=dim(Rinf)[1], ncol=1))
  tem1 = 1 - colSums(Rinf^2)
  w<-backsolve(U,one,transpose=T)
  FRr<-crossprod(w,Rinf)
  FRF<-crossprod(w)
  tem2<-(1-FRr)^2/as.numeric(FRF)
  results<-(tem1+tem2)*((n-1)*sigma+0)/(n-1+0)
  rett<-matrix(results,ncol=1,byrow=F)
  return(rett)
}
