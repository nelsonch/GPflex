#' Compute Four components of a GP with a matern structure
#'
#' The function computes the posterior probability at log scale without the normalization constant, the correlation matrix of 
#' points in the training set, the estiamted constant mean and the variance of the GP. 
#' 
#' @param y True outputs of x. Size is n. 
#' @param lambda sensitivity parameters at lambda scale.
#' @param d a n by n matrix with all elements 0. Must be passed from dist.R.matern(). 
#' @param odut the indices by which to place the n*(n-1)/2 correlation 
#' in the correlation matrix filled by rows. Must be passed from dist.R.matern().  
#' @param dup a n*(n-1)/2 by d matrix contains the n*(n-1)/2 combinations 
#' of (x_ik - x_jk) for any i < j and k < d. Must be passed from dist.R.matern().  
#' @param cp smoothness parameters one of 0, 1, 2, >=3. 
#' @return R returns the correlation matrix. 
#' @return logpost returns the posterior probability at log scale without the normalization constant.  
#' @return mu returns the estimated constant mean of the GP. 
#' @return sigma returns the estimated variance of the GP.  


log.post.matern<-function(d, odut,dup, y, lambda, cp){
  ##rho must be a  vector
  #####distance matrix & R  
  y=y; 
  rho=exp(-lambda)/(1+exp(-lambda))
  phi=-4*log(rho); 
  #phi=-4*log(rho)
  d<-d; odut<-odut; dup<-dup
  n<-dim(d)[1]
  ddp<-length(cp)
  dis<-rep(1, dim(dup)[1])
  for(j in 1:ddp){
    h <- phi[j]*dup[,j] 
    tem<-cp[j]
    if(tem==0){
      dis <- dis * exp(-h)
    }
    
    if(tem==1){
      dis <- dis *exp(-h)*(h+1)
    }
    
    if(tem==2){
      dis <- dis *exp(-h) * ((h^2/3)+h+1)  
    }
    
    if(tem>=3){
      dis <- dis * exp(-phi[j]*(dup[,j])^2) 
    }
  }
  dis.m<-as.matrix(dis)
  d[odut]<-dis.m
  d<-d+t(d)
  diags = seq(0,n*(n-1),by=n) + 1:n
  d[diags] = 1
  R<-d
  U<-try(chol(R),silent=T)
  if (class(U)=="try-error")
  { 
    logpost=-9e99; mu.hat=logpost; sigma.hat=logpost
  }
  else{
    F1=mat.or.vec(n,1)+1
    k=1
    S=backsolve(U,F1,transpose=T)
    S1=crossprod(S)
    G=backsolve(U,y,transpose=T)
    A=crossprod(S,G)
    ##mu.hat
    mu.hat=solve(S1,A)
    B=y-mu.hat*F1
    B1=backsolve(U,B,transpose=T)
    B2=crossprod(B1)
    sigma.hat=(1/(n-k))*B2
    logdet = sum(log(diag(U)))
    S2=log(det(crossprod(S)))
    tem1=0.5*(exp(-lambda)/(1+exp(-lambda)))^(-0.5)*(1+exp(-lambda))^(-2)*exp(-lambda)
    tem11<-sum(log(tem1))
    logprior<-tem11
    logpost=logprior-(0+(n-k)/2)*log(0+(n-k)*sigma.hat/2)-0.5*S2-logdet
  }
  para=list(logpost=logpost,R=R,mu=mu.hat,sigma=sigma.hat)
  return(para)
}
