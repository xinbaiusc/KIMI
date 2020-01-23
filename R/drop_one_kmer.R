#' @title Drop one kmer.
#' @description Select the kmer having the least likelihood of being relevant by WMW test.
#' @param kmerfX_normal The normalized kmer frequency matrix.
#' @param y The label vector, taking values of 0 or 1 for two classes of sequences.
#' @return drop_idx The index of kmer to be dropped.
#' @export

drop_one_kmer <- function(kmerfX_normal, y)
{

y0_idx <- which(y == 0)
y1_idx <- which(y == 1)
p <- dim(kmerfX_normal)[2]
pvalue=rep(0,p)
for (col in 1:p)
{
pvalue[col] <- wilcox.test(kmerfX_normal[y0_idx, col], kmerfX_normal[y1_idx, col],paired = FALSE)$p.value
}
drop_idx <- which(pvalue==max(pvalue))[1]
drop_idx

}