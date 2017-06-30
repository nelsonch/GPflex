#'The function is the R equivalent of repmat (matlab)
#' 
#' @param X is a mtrix. 
#' @param m is the number of repeats in the direction of row. 
#' @param n is the number of repeats in the direction of column. 
#' @return the resulting matrix. 
repmat = function(X, m, n){
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx*n)), mx*m, nx*n, byrow=T)
}