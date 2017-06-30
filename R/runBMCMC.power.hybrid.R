#' Gaussian process with a hybrid power-exponential (Power-Hybrid) correlation structure.
#'
#' The power-exponential structure is
#' \deqn{ R(x_i, x_j)= \exp(\sum_{k=1}^{d} \theta_k (x_{i, k} - x_{j,k})^{\alpha_k}),}
#' where \eqn{ 1 \le  i, j \le n}, \eqn{\theta_k > 0 } is the sensitivity parameter and
#' \eqn{ 1 \le \alpha_k \le 2 } is the smoothness parameter.
#' \deqn{ \lambda_k = \log_{e}(\rho_k/(1-\rho_k)),}
#' where \eqn{\rho_k = \exp(-\theta_k/4)}. After the above \eqn{\lambda} transformation, the sensitivity parameter
#' is unconstrained.
#' @param nmcmc Length of the MCMC chain. Default is 10,000.
#' @param burn Burn-in period. Default is 4,000.
#' @param thin Thinning parameter. Default is 10.
#' @param x Training set. Either a matrix or a data frame.
#' @param y True outputs of x.
#' @param xtest Hold-out set, of which the ouputs will be predicted.
#' @param lambda.ini Initial values of the sensitivity parameters of the Power-Hybrid structure, at \eqn{\lambda} scale.
#' @param lambda.w.ini Initial values of the lengths (variances) for adaptive MCMC in estimating the sensitivity parameters of the Power-Hybrid structure, at \eqn{\lambda} scale.
#' @param p Smoothness parameters from 1 to 2.
#' @return pred.y return the predicted mean.
#' @return pred.var return the predicted variance.
#' @return accept.rate.lambda return a vector (1*d) of the acceptance rate of lambda.
#' @return mcmc.matrix.lambda return a matrix (nmcmc * d) of the sampled lambda.
#'
runBMCMC.power.hybrid<-function(nmcmc = 10000, burn = 4000,
                                thin = 10, x, y, xtest, lambda.ini,
                                lambda.w.ini, p){

  dp <- ncol(x)
  ##########
  mcmc.ma.lambda <- matrix(nrow=nmcmc, ncol=dp, byrow=T)
  accept.lambda <- matrix(c(rep(0, dp)), nrow=1, ncol=dp, byrow=T)
  res <- matrix(nrow=(nmcmc-burn)/thin, ncol=nrow(xtest))
  v.term2 <- matrix(nrow=(nmcmc-burn)/thin, ncol=nrow(xtest))
  ###########
  j <- 0
  ff <- dist.R.power.hybrid(x, p)
  cov.dd <- cov.r1.dis.power.hybrid(x, xtest, p)

  lambda = lambda.ini
  for(i in 1 : nmcmc){
    ####### The first two iterations. variance is not computed.
    ####### We use lambda.w.ini instead.
    if(i <= 2){
       for(k in 1:dp){
        lambda.cond <- NULL
        lambda.cond <- lambda
        temlambda <- 0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
        lambda.cond[k] <- temlambda
        com.phi<-log.post.power.hybrid(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond)$logpost-
        log.post.power.hybrid(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda)$logpost
        u<-runif(1)
        if(log(u) < com.phi){
        lambda[k] <- lambda.cond[k]
        if(i > burn){
          accept.lambda[1,k]=accept.lambda[1,k]+1
        }
      }
    }
    mcmc.ma.lambda[i,] <- lambda

    if( i > burn && ((i-burn)%%thin==0)){
      j = j + 1
      fit10 <- pred1.power.hybrid(ff$d, ff$odut, ff$dup, cov.dd$dis, y,lambda)
      res[j,] <- fit10$res
      v.term2[j,] <- fit10$v.term2
    }

    if ((i %% (0.1*nmcmc))==0){
      print(c(i/nmcmc))
    }
    }else{
    for(k in 1:dp){
      lambda.cond <- NULL
      lambda.cond <- lambda
      tg <- mcmc.ma.lambda[1:(i-1), k]
      temlambda <- 0.95*rnorm(1, mean=lambda[k], sd=2.38*sd(tg))+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k] <- temlambda
      com.phi <- log.post.power.hybrid(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond)$logpost-
        log.post.power.hybrid(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda)$logpost
      u <- runif(1)
      if(log(u) < com.phi){
        lambda[k] <- lambda.cond[k]
        if(i > burn){
        accept.lambda[1,k]=accept.lambda[1,k]+1
        }
      }
    }
    mcmc.ma.lambda[i,] <- lambda

    if( i > burn && ((i-burn)%%thin==0)){
      j = j + 1
      fit10 <- pred1.power.hybrid(ff$d,ff$odut,ff$dup,cov.dd$dis, y,lambda)
      res[j,] <- fit10$res
      v.term2[j,] <- fit10$v.term2
    }

    if ((i%%(0.1*nmcmc))==0){
      print(c(i/nmcmc))
    }
    }
  }
  var.term1 <- apply(res, 2, var)
  var.term2 <- apply(v.term2, 2, mean)
  pred.var <- var.term1+var.term2
  pred.y <- apply(res, 2, mean)
  accept_rate <- accept.lambda/(nmcmc-burn)
  m<-list(pred.y = pred.y, pred.var = pred.var, accept.rate.lambda = accept_rate,
          mcmc.matrix.lambda = mcmc.ma.lambda)
  return(m)
}
