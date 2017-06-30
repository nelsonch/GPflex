#' Gaussian process with a hybrid matern (Matern-Hybrid) correlation structure.
#'
#' @param nmcmc Length of the MCMC chain. Default is 10,000.
#' @param burn Burn-in period. Default is 4,000.
#' @param thin Thinning parameter. Default is 10.
#' @param x Training set. Either a matrix or a data frame.
#' @param y True outputs of x.
#' @param xtest Hold-out set, of which the ouputs will be predicted. Must have the same dimension as x.
#' @param lambda.ini Initial values of the sensitivity parameters, at \eqn{\lambda} scale.
#' @param lambda.w.ini Initial values of the lengths for adaptive MCMC in estimating the sensitivity parameters, at \eqn{\lambda} scale.
#' @param cp Smoothness parameters one of 0, 1, 2, >=3. Must be a vector with the same dimension as x.
#' @return pred.y return the predicted mean.
#' @return pred.var return the predicted variance.
#' @return accept.rate.lambda return a vector (1*d) of the acceptance rate of lambda.
#' @return mcmc.matrix.lambda return a matrix (nmcmc * d) of the sampled lambda.
runBMCMC.matern.hybrid<-function(nmcmc = 10000,
                                 burn = 4000,
                                 thin = 10,
                                 x, y, xtest, lambda.ini, lambda.w.ini, cp){

  mcmc.ma.lambda<-NULL
  mcmc.ma.lambda<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.lambda<-NULL
  accept.lambda<-matrix(c(rep(0,dp)),nrow=1,ncol=dp,byrow=T)
  lambda<-lambda.ini
  lambda.w.ini<-lambda.w.ini

  res<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest))
  v.term2<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest))
  ff<-dist.R.matern(x)
  cov.dd<-cov.r1.dis.matern(x,xtest)
  j <- 0

  lambda = lambda.ini
  for(i in 1 : nmcmc){
    ####### The first two iterations. variance is not computed.
    ####### We use lambda.w.ini instead.
    if(i <= 2){
      for(k in 1:dp){
        lambda.cond<-NULL
        lambda.cond<-lambda
        temlambda<-0.95*rnorm(1, mean=lambda[k], sd = 2.38 * lambda.w.ini[k])+
          0.05*rnorm(1, mean=lambda[k], sd = 0.1)
        lambda.cond[k]<-temlambda
        com.phi<-log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond, cp)$logpost-
          log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda, cp)$logpost
        u<-runif(1)
        if(log(u)<com.phi){
          lambda[k]<-lambda.cond[k]
          accept.lambda[1,k]=accept.lambda[1,k]+1
        }
      }
      mcmc.ma.lambda[i,]<-lambda

      if(i>burn&&((i-burn)%%thin==0)){
        j=j+1
        fit10<-pred1.matern(ff$d,ff$odut,ff$dup,cov.dd$dis, y,lambda, cp)
        res[j,]=fit10$res
        v.term2[j,]=fit10$v.term2
      }
      if ((i%%(0.1*nmcmc))==0){
        print(c(i/nmcmc))
      }
    }else{

      ###### lambda
      for(k in 1:dp){
        lambda.cond<-NULL
        lambda.cond<-lambda
        tg<-mcmc.ma.lambda[1:(i-1), k]
        temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*sd(tg))+
          0.05*rnorm(1, mean=lambda[k], sd=0.1)
        lambda.cond[k]<-temlambda
        com.phi<-log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond, cp)$logpost-
          log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda, cp)$logpost
        u<-runif(1)
        if(log(u)<com.phi){
          lambda[k]<-lambda.cond[k]
          accept.lambda[1,k]=accept.lambda[1,k]+1
        }
      }

      mcmc.ma.lambda[i,]<-lambda

      if(i>burn&&((i-burn)%%thin==0)){
        j=j+1
        fit10<-pred1.matern(ff$d,ff$odut,ff$dup,cov.dd$dis, y,lambda, cp)
        res[j,]=fit10$res
        v.term2[j,]=fit10$v.term2
      }

      if ((i%%(0.1*nmcmc))==0){
        print(c(i/nmcmc))
      }

    }
  }

  var.term1<-apply(res,2,var)
  var.term2<-apply(v.term2,2,mean)
  pred.var<-var.term1+var.term2
  pred.y<-apply(res,2,mean)
  accept_rate <- accept.lambda/(nmcmc-burn)

  m<-list(pred.y=pred.y, pred.var=pred.var, accept.rate.lambda=accept_rate,
          mcmc.matrix.lambda=mcmc.ma.lambda)
  return(m)
}
