#' Compute part of the correlation matrix based on a matern structure
#' 
#' For every point in the test set, the function computes the distances of it to all points in the training set, i.e.,
#' (x_hk - x_ik) for h <= m (size of the test set), j <=n (size of the training set), and k < d.
#'
#' @param x Training set. Either a matrix or a data frame. The same size is n and dimension is d.
#' @param xtest Hold-out set, of which the ouputs will be predicted.
#' @return dis.tem returns a m*n by d matrix contains the m*n combinations of (x_hk - x_jk) for h <= m, j <=n, and k < d.

cov.r1.dis.matern<-function(x, xtest){
  n<-dim(x)[1]
  m<-dim(xtest)[1]
  inds=n*m
  indi=repmat(matrix(c(1:n),ncol=n),1,m)
  indj=repmat(matrix(c(1:m),ncol=m),n,1)
  indj=t(unlist(list(indj)))
  dis.tem<-abs(x[indi,]-xtest[indj,])
  return(list(dis=dis.tem))
}