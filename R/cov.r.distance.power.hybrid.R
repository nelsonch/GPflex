#' Compute part of the correlation matrix based on a PowExp structure
#' 
#' For every point in the test set, the function computes the powered distances of it to all points in the training set
#' without multiplied by the sensitivity parameters, i.e.,
#' (x_hk - x_ik)^p_k for h <= m (size of the test set), j <=n (size of the training set), and k < d.
#'
#' @param x Training set. Either a matrix or a data frame. The same size is n and dimension is d.
#' @param xtest Hold-out set, of which the ouputs will be predicted.
#' @param p Smoothness parameters from 1 to 2.
#' @return dis returns a m*n by d matrix contains the m*n combinations of (x_hk - x_jk)^p_k for h <= m, j <=n, and k < d.

cov.r1.dis.power.hybrid <- function(x, xtest, p){
  n<-dim(x)[1]
  m<-dim(xtest)[1]
  inds=n*m
  indi = repmat(matrix(c(1:n),ncol=n),1,m)
  indj = repmat(matrix(c(1:m),ncol=m),n,1)
  indj=t(unlist(list(indj)))
  dis.tem<-(x[indi,]-xtest[indj,])
  dis<-NULL
  for(i in 1:dim(dis.tem)[2]){
    dis<-cbind(dis, (abs((dis.tem[,i])))^(p[i]))
  }
  return(list(dis=dis))
}
