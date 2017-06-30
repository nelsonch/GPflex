#' Gaussian process with a full power-exponential (Full-Hybrid) correlation structure.
#'
#' The power-exponential structure is
#' \deqn{ R(x_i, x_j)= \exp(\sum_{k=1}^{d} \theta_k (x_{i, k} - x_{j,k})^{\alpha_k}),}
#' where \eqn{ 1 \le  i, j \le n}, \eqn{\theta_k > 0 } is the sensitivity parameter and
#' \eqn{ 1 \le \alpha_k \le 2 } is the smoothness parameter.
#' \deqn{ \lambda_k = \log_{e}(\rho_k/(1-\rho_k)),}
#' where \eqn{\rho_k = \exp(-\theta_k/4)}. After the above \eqn{\lambda} transformation, the sensitivity parameter
#' is unconstrained.
#' \deqn{ \gamma_k = -\exp(\frac{2-\alpha_k}{\alpha_k-1}).}
#' The smoothness parameter is also unconstrained at \eqn{\gamma} scale.
#'
#' @param nmcmc Length of the MCMC chain. Default is 10,000.
#' @param burn Burn-in period. Default is 4,000.
#' @param thin Thinning parameter. Default is 10.
#' @param x Training set. Either a matrix or a data frame.
#' @param y True outputs of x.
#' @param xtest Hold-out set, of which the ouputs will be predicted. Must have the same dimension as x.
#' @param lambda.ini Initial values of the sensitivity parameters, at \eqn{\lambda} scale.
#' @param lambda.w.ini Initial values of the lengths (variances) for adaptive MCMC in estimating the sensitivity parameters, at \eqn{\lambda} scale.
#' @param gamma.ini Initial values of the smoothness parameters, at \eqn{\gamma} scale.
#' @param gamma.w.ini Initial values of the lengths (variances) for adaptive MCMC in estimating the smoothness parameters, at \eqn{\gamma} scale.
#' @return pred.y return the predicted mean.
#' @return pred.var return the predicted variance.
#' @return accept.rate.lambda return a vector (1*d) of the acceptance rate of lambda.
#' @return accept.rate.gamma return a vector (1*d) of the acceptance rate of gamma.
#' @return mcmc.matrix.lambda return a matrix (nmcmc * d) of the sampled lambda.
#' @return mcmc.matrix.gamma return a matrix (nmcmc * d) of the sampled gamma.
#'
runBMCMC.power.full<-function(nmcmc=10000, burn=4000, thin=10, x, y, xtest, lambda.ini,
                              lambda.w.ini, gamma.ini, gamma.w.ini){

  dp<-ncol(x)
  j=0
  ###### lambda
  mcmc.ma.lambda<-NULL
  mcmc.ma.lambda<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.lambda<-NULL
  accept.lambda<-matrix(c(rep(0,dp)),nrow=1,ncol=dp,byrow=T)
  ###### gamma
  mcmc.ma.gamma<-NULL
  mcmc.ma.gamma<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.gamma<-NULL
  accept.gamma<-matrix(c(rep(0,dp)),nrow=1,ncol=dp,byrow=T)

  lambda <- lambda.ini
  gamma <- gamma.ini

  res<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest))
  v.term2<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest))

  ff<-dist.R.power.full(x)
  cov.dd<-cov.r1.dis.power.full(x,xtest)

  for(i in 1 : nmcmc){
    ####### The first two iterations. variance is not computed.
    ####### We use lambda.w.ini and gamma.w.ini instead.
    if(i <= 2){
      for(k in 1:dp){
        lambda.cond<-NULL
        lambda.cond<-lambda
        temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
          0.05*rnorm(1, mean=lambda[k], sd=0.1)
        lambda.cond[k]<-temlambda
        com.phi<-log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y, lambda.cond, gamma)$logpost-
          log.post.power.full(d=ff$d, odut=ff$odut, tem.dup=ff$tem.dup, y, lambda, gamma)$logpost
        u<-runif(1)
        if(log(u)<com.phi){
          lambda[k]<-lambda.cond[k]
          accept.lambda[1,k]=accept.lambda[1,k]+1
        }

        ###########gamma
        gamma.cond<-NULL
        gamma.cond<-gamma
        temgamma<-0.95*rnorm(1, mean=gamma[k], sd=2.38*gamma.w.ini[k])+
          0.05*rnorm(1, mean=gamma[k], sd=0.1)
        temg<-(temgamma)
        if((temg >=-20) & (temg <=20)){
          gamma.cond[k]<-temg
          com.gamma<-log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma.cond)$logpost-
            log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
          u<-runif(1)
          if(log(u)<com.gamma){
            gamma[k]<-gamma.cond[k]
            accept.gamma[1,k]=accept.gamma[1,k]+1
          }
        }
      }
      mcmc.ma.lambda[i,] <- lambda
      mcmc.ma.gamma[i,] <- gamma
      if(i>burn&&((i-burn)%%thin==0)){
        j=j+1
        fit10<-pred1.power.full(ff$d,ff$odut,ff$tem.dup,cov.dd$dis.tem, y, lambda, gamma)
        res[j,]=fit10$res
        v.term2[j,]=fit10$v.term2
      }
      if ((i%%(0.1*nmcmc))==0){
        print(c(i/nmcmc))
      }
    }else{
      for(k in 1:dp){
        lambda.cond<-NULL
        lambda.cond<-lambda
        tg<-mcmc.ma.lambda[1:(i-1), k]
        temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*sd(tg))+
          0.05*rnorm(1, mean=lambda[k], sd=0.1)
        lambda.cond[k]<-temlambda
        com.phi<-log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda.cond, gamma)$logpost-
          log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
        u<-runif(1)
        if(log(u)<com.phi){
          lambda[k]<-lambda.cond[k]
          accept.lambda[1,k]=accept.lambda[1,k]+1
        }

        #########gamma
        gamma.cond<-NULL
        gamma.cond<-gamma
        tgg<-mcmc.ma.gamma[1:(i-1), k]
        temgamma<-0.95*rnorm(1, mean=gamma[k], sd=2.38*sd(tgg))+
          0.05*rnorm(1, mean=gamma[k], sd=0.1)
        temg<-(temgamma)
        if((temg >=-20) & (temg <=20)){
          gamma.cond[k]<-temg
          com.gamma<-log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup, y,lambda, gamma.cond)$logpost-
            log.post.power.full(d=ff$d, odut=ff$odut,tem.dup=ff$tem.dup,y,lambda, gamma)$logpost
          u<-runif(1)
          if(log(u)<com.gamma){
            gamma[k]<-gamma.cond[k]
            accept.gamma[1,k]=accept.gamma[1,k]+1
          }
        }
      }
      mcmc.ma.lambda[i,]<-lambda
      mcmc.ma.gamma[i,]<-gamma
      if(i>burn&&((i-burn)%%thin==0)){
        j=j+1
        fit10<-pred1.power.full(ff$d,ff$odut,ff$tem.dup,cov.dd$dis.tem, y, lambda, gamma)
        res[j,]=fit10$res
        v.term2[j,]=fit10$v.term2
      }
      if ((i%%(0.1*nmcmc))==0){
        print(c(i/nmcmc))
      }
    }
  }

  ##########
  var.term1<-apply(res, 2, var)
  var.term2<-apply(v.term2, 2, mean)
  pred.var<-var.term1 + var.term2
  pred.y<-apply(res, 2, mean)

  #########acceptance rate
  accept_rate_lambda <- accept.lambda/(nmcmc-burn)
  accept_rate_gamma <- accept.gamma/(nmcmc-burn)

  m<-list(pred.y=pred.y, pred.var=pred.var, accept.rate.lambda=accept_rate_lambda,
          mcmc.matrix.lambda=mcmc.ma.lambda,
          accept.rate.gamma=accept_rate_gamma,
          mcmc.ma.gamma=mcmc.ma.gamma)
  return(m)
}


