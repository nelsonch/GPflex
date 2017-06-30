#' Compute part of the correlation matrix based on a PowExp structure
#'
#' The function computes the distances of any two points in the training set, i.e.,  
#' (x_ik - x_jk) for any i < j and k < d. 
#' 
#' @param x Training set. Either a matrix or a data frame. The same size is n and dimension is d.  
#' @return d returns a n by n matrix with all elements 0. 
#' @return odut returns the indices by which to place the n*(n-1)/2 correlation in the correlation matrix filled by rows. 
#' @return tem.dup returns a n*(n-1)/2 by d matrix contains the n*(n-1)/2 combinations of (x_ik - x_jk) for any i < j and k < d.   

dist.R.power.full<-function(x){
  n<-nrow(x)
  d<-mat.or.vec(n,n)
  inds=n*(n-1)/2
  indi=mat.or.vec(inds,1)
  indj=mat.or.vec(inds,1)
  ind=1
  for (ii in 2:n-1) 
  {
    indi[ind:(ind+n-ii-1)]=ii 
    indj[ind:(ind+n-ii-1)]=(ii+1):n
    ind=ind+n-ii 
  }
  tem.dup<-(x[indi,]-x[indj,])
  odut=indi + n*(indj-1)
  ttm<-list(d=d, odut=odut, tem.dup=tem.dup)
  return(ttm)
}