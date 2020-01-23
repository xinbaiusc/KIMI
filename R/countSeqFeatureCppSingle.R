#' @title Count kmer from sequences from single direction.
#' @description Count kmer from sequences from single direction using the cpp file.
#' @param RseqDNA The DNA sequence.
#' @param k The kmer size.



countSeqFeatureCppSingle <- function(RseqDNA, k) {
    .Call('KIMI_countSeqFeatureCppSingle', PACKAGE = 'KIMI', RseqDNA, k)
}