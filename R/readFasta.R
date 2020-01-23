#' @title Read dasta.
#' @description Reads a fasta/fastq format file containing contigs, count kmer frequency for each contig, and then output a kmer frequency matrix.
#' @param FastaDir The input directory of the fasta/fastq file.
#' @param w The kmer size.
#' @param is_single_strand A bool variable indicating counting single strand or double strand, default to be TRUE.
#' @return A kmer frequency matrix. Each row corresponds to a contig. Each column corresponds to a kmer.
#' @export


readFasta <- function(FastaDir, w, is_single_strand = TRUE)
{

lineNum <- 0
con <- file(FastaDir, open = "r")

while ( length(line <- readLines(con, n = 1)) > 0 ) 
{
if( toupper(strsplit(line, "")[[1]][1]) == "A" | toupper(strsplit(line, "")[[1]][1]) == "C" | toupper(strsplit(line, "")[[1]][1]) == "G" | toupper(strsplit(line, "")[[1]][1]) == "T" )
{
lineNum <- lineNum + 1
}
}

n <- lineNum

close(con)

if (is_single_strand)
{
kmerfX <- matrix(0, n, 4^w)
}else if( w %% 2 == 0 ){
kmerfX <- matrix(0, n, (4^w+2^w)/2)
}else{
kmerfX <- matrix(0, n, (4^w)/2)
}

lineNum <- 0
con <- file(FastaDir, open = "r")
while ( length(line <- readLines(con, n = 1)) > 0 ) 
{
if( toupper(strsplit(line, "")[[1]][1]) == "A" | toupper(strsplit(line, "")[[1]][1]) == "C" | toupper(strsplit(line, "")[[1]][1]) == "G" | toupper(strsplit(line, "")[[1]][1]) == "T" )
{
lineNum <- lineNum + 1
seq_list=strsplit(line, "")
seq_vector = as.vector(unlist(seq_list))

if (is_single_strand){
kmerfX[lineNum, ] <- as.matrix(countSeqFeatureCppSingle(seq_vector, w)$kmerCount)
}else{
kmerfX[lineNum, ] <- as.matrix(countSeqFeatureCpp(seq_vector, w)$kmerCount)
}
}
}  
close(con)

kmerfX
}