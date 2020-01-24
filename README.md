# KIMI: R package based on Knockoff Inference for Motif Identification from molecular sequences with false discovery rate control 

Version: 1.0.0

Authors: Xin Bai, Jie Ren, Yingying Fan, Fengzhu Sun

Maintainer: Xin Bai, <xinbai@usc.edu>

## Description

The R package provides functions based on knockoff inference to identify motifs from binary types of molecular sequences stored in a fasta file, such as the assembled contigs from metagenomic data. 

For each query sequence, KIMI first counts  the number of occurrences of all k-mers by a c++ program using a hash table. The kmer counts are then normalized into frequencies. For each particular kmer, 
the kmer frequency vector of all sequences are further normalized with mean 0 and standard diviation 1. Then the WMW test is applied for testing the significance of each kmer being relevant, and
the kmer that is least likely to be relevant is dropped. We then apply knockoff inference for selecting relevant kmers . 

## Dependencies

R packages "VirFinder", "glmnet" and "Rcpp" are needed to be installed before Installation of KIMI.

To install "glmnet" and "Rcpp", start R and enter,

<pre>
install.packages("glmnet", dependencies=TRUE)
install.packages("Rcpp", dependencies=TRUE)
</pre>

To install "VirFinder", please refer to 
[VirFinder](https://github.com/jessieren/VirFinder)

## Installation

To install the R package KIMI, follow the instuctions on

<http://cran.r-project.org/doc/manuals/r-release/R-admin.html#Installing-packages>

## Usage

The KIMI R package takes a fasta file containing binary types of sequences, sequence labels, and the prespecified kmer size as inputs, and outputs indexes of selected kmers sorted according to the alphabet order.

For example, one can first read the sequences and count the occurrence of 5-mers,

<pre>
library(KIMI)
F <- readFasta(&lt;path_to_the_fasta_file&gt;, w = 5)
</pre>
