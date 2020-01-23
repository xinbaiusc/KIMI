#' @title Normalize into standard gaussian.
#' @description Normalize each column of the kmer frequency matrix into a vector with mean 0 and standard deviation 1.
#' @param kmerfX The kmer frequency matrix.
#' @return A normalized kmer frequency matrix. Each column has mean 0 and standard deviation 1.
#' @export
 
normalize_standard_gaussian <- function(kmerfX)
{

kmerfX_normal <- kmerfX
for (col in 1:dim(kmerfX)[2])
{
if (sum(kmerfX[,col]^2) > 0){
kmerfX_normal[,col] <- kmerfX_normal[,col]/sqrt(sum(kmerfX_normal[,col]^2))
}
}

kmerfX_normal
}