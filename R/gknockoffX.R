#' @title Design matrix generation.
#' @description Function generating the knockoff design matrix.
#' @param X The normalized kmer frequency matrix with one column dropped.
#' @param Omega The estimated precision matrix.
#' @return n x (2p) design matrix = [X X_ko]. X_ko is the knockoff matrix.
#' @export

gknockoffX <- function(X, Omega){
  n = nrow(X)
  p = ncol(X)
  
  # calculate minimum eigenvalue of Sigma = inv(Omega), i.e. 1/max eigenvalue of Omega
  r = eigen(Omega)
  max.eig = r$values[1]
 
  # diagonal matrix for construction of knockoff variables
  s = (1/max.eig)*diag(p)
  obj = sqrtm(2*s - s%*%Omega%*%s)
  B = obj$B
  A = diag(p) - s%*%Omega
  # now construct knockoff variables conditional on X
  X_ko = X%*%A + matrix(rnorm(n*p),n,p)%*%B
  Xnew = cbind(X, X_ko)
  
  return(Xnew)
}