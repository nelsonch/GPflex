bceMCMC.matern<-function(nmcmc,burn,thin,x,y,xtest1, lambda.ini, lambda.w.ini, cp){
  nmcmc=nmcmc
  burn=burn
  thin=thin
  x=x
  xtest1=xtest1
  y=y
  dp<-ncol(x)
  j=0
  mcmc.ma.lambda<-NULL
  mcmc.ma.lambda<-matrix(nrow=nmcmc,ncol=dp,byrow=T)
  accept.lambda<-NULL
  accept.lambda<-matrix(c(rep(0,dp)),nrow=1,ncol=dp,byrow=T)
  reasonable.lambda<-NULL
  reasonable.lambda<-matrix(c(rep(0,dp)),nrow=1,ncol=dp,byrow=T)
  lambda<-lambda.ini
  lambda.w.ini<-lambda.w.ini
  res<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest1))
  v.term2<-matrix(nrow=(nmcmc-burn)/thin,ncol=nrow(xtest1))
  ff<-dist.R.matern(x)
  cov.dd<-cov.r1.dis.matern(x,xtest1)
  for(i in 1:1){
    for(k in 1:dp){
      lambda.cond<-NULL
      lambda.cond<-lambda  
      temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k]<-temlambda
      com.phi<-log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond, cp)$logpost-
        log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda, cp)$logpost
      u<-runif(1)
      if(log(u)<com.phi){
        lambda[k]<-lambda.cond[k]
        accept.lambda[1,k]=accept.lambda[1,k]+1
      }
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
    }
    mcmc.ma.lambda[i,]<-lambda
  }
  
  for(i in 2:2){
    for(k in 1:dp){
      lambda.cond<-NULL
      lambda.cond<-lambda  
      temlambda<-0.95*rnorm(1, mean=lambda[k], sd=2.38*lambda.w.ini[k])+
        0.05*rnorm(1, mean=lambda[k], sd=0.1)
      lambda.cond[k]<-temlambda
      com.phi<-log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup, y,lambda.cond, cp)$logpost-
        log.post.matern(d=ff$d, odut=ff$odut,dup=ff$dup,y,lambda, cp)$logpost
      u<-runif(1)
      if(log(u)<com.phi){
        lambda[k]<-lambda.cond[k]
        accept.lambda[1,k]=accept.lambda[1,k]+1
      }
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
    }
    mcmc.ma.lambda[i,]<-lambda
  }
  
  for(i in 3:nmcmc){
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
      reasonable.lambda[1, k]=reasonable.lambda[1,k]+1
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
  m<-list(reasonable.lambda=reasonable.lambda, accept.lambda=accept.lambda, 
          mcmc.ma.lambda=mcmc.ma.lambda, res=res, v.term2=v.term2)
  return(m)
}