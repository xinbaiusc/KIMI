#' @title KIMI for motif identification.
#' @description KIMI for motif identification based on the knockoff plus statistic.
#' @param X The normalized kmer frequency matrix with one column dropped.
#' @param y The label vector, taking values of 0 or 1 for two classes of sequences.
#' @param q The target FDR.
#' @param alpha The parameter for the weight of L1 penalty.
#' @return The set of indexes of selected kmers.
#' @export

KIMI_plus <- function(X, y, q, alpha = 1){
  
  n = nrow(X)
  p = ncol(X)
  
  # calculate the precision matrix 
  if (n < p)
  {
  stop(" n must be greater than p in order to make the covariance matrix invertible! ")
  }else{
  Y = solve(cov(X))
  Omega = Y/abs(median(Y))
  X = X/sqrt(abs(median(Y)))
  }
  
  # calculate the knockoff variables
  Xnew = gknockoffX(X, Omega)
  
  # calculate the knockoff statistics using glmnet
  W = stat.glmnet_coefdiff(X, Xnew[,(p+1):(2*p)], y, nfolds = 10, family = "binomial", alpha)
  
  # find the knockoff threshold T
  t = sort(c(0, abs(W)))
  ratio = c(rep(0, p))
  for (j in 1:p) {
    ratio[j] = (1+sum(W <= -t[j]))/max(1, sum(W >= t[j]))
  }
  id = which(ratio <= q)[1]
  if(length(id) == 0){
    T = Inf
  } else {
    T = t[id]
  }
  
  # set of discovered variables
  S_plus = which(W >= T)
  
  return(S_plus)
}