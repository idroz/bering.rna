% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rna.mattest.R
\name{rna.mattest}
\alias{rna.mattest}
\title{Array-based t-test}
\usage{
rna.mattest(expr, ref, samp, p.value = NULL)
}
\arguments{
\item{expr}{numeric matrix or data.frame consisting of gene expression values}

\item{ref}{numeric vector. Index of reference samples in the array}

\item{samp}{numeric vector. Index of case samples in the array}

\item{p.value}{double or NULL (default). If NULL returns a vector of p-values
for each genes, otherwise returns a logical vecotr of genes with
p-values below the specified threshold.}
}
\description{
Perform two-sample t-test to evaluate differential expression of genes from
two experimental conditions or phenotypes
}

