#' Compute the correlation matrix based on a matern structure
#' 
#' The function computes the m by n correlation between points in the hold set and points in the training set. 
#' 
#' @param dis.tem a m*n by d matrix contains the m*n combinations of (x_hk - x_jk) for h <= m, j <=n, and k < d. 
#' Must be passed from cov.r1.dis.matern().   
#' @param lambda sensitivity parameters at lambda scale.
#' @param n size of training set. 
#' @param m size of test set. 
#' @param cp smoothness parameters one of 0, 1, 2, >=3.
#' @return r.matrix returns the correlation matrix between all points in the test set and all points in the training set.
#'  
cov.r1.matern<-function(dis,lambda, n, m, cp){
  rho=exp(-lambda)/(1+exp(-lambda))
  phi=-4*log(rho); 
  dis<-dis
  dis.t<-rep(1, dim(dis)[1])
  ddp<-length(cp)
  for(j in 1:ddp){
    tem<-cp[j]
    h <- phi[j]*dis[,j]
    
    if(tem==0){
      dis.t <- dis.t * exp(-h)
    }
    
    if(tem==1){
      dis.t <- dis.t *exp(-h) * (h+1)
    }
    
    if(tem==2){
      dis.t <- dis.t * exp(-h)*(((h)^2/3)+h+1)  
    }
    
    if(tem>=3){
      dis.t <- dis.t * exp(-phi[j]*(dis[,j])^2) 
    }
  }
  r<-dis.t
  r.matrix<-matrix(r,nrow=m,ncol=n,byrow=T)
  return(r.matrix)
}
