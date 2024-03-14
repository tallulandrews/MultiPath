\name{get_pseudobulk}
\alias{get_pseudobulk}
\title{Calculate Pseudobulks}
\description{
	Calculates pseudobulks given a count matrix and cell-type and sample labels.
}
\usage{
	get_pseudobulks(mat, clusters, donors, method=c("sum", "mean"), trim=10, refactor=TRUE)
}
\arguments{
  \item{mat}{a matrix of umi counts or normalized expression for each cell, supports sparse matrices.}
  \item{clusters}{a vector of cluster or cell-type labels for each cell.}
  \item{donors}{a vector of donor or sample labels for each cell.}
  \item{method}{whether to add up or average the expression in each cell-type x sample group of cells.}
  \item{trim}{cell-type x sample groups with fewer than this many cells will not be included in the pseudobulk matrix}
  \item{refactor}{whether to factor the clusters & donors vectors - recommended if trimming}
}
\details{ 
	Calculates pseudobulk expression for performing differential expression (DE). If using a negative-binomial model based DE method such as edgeR the "sum" of the raw umi counts should be used. If using a Gaussian or non-parametric DE method then the "mean" of normalized expression should be calculated.

	To avoid noise from samples where only a small number of cells of a particular cell-type are found, it automatically filters out cases where there are fewer than `trim` cells of a particular cell-type in a particular sample.
	
	Note: to avoid issues retreiving sample & cell-type IDs from the output column names, all underscores in the original clusters and donor ids are replaced with a dash("-").
}
\value{
	A matrix of pseudobulk expression with column names in the format: [cluster]_[donor].
}
\examples{
  example_data <- generate_test_cellcounts()
  test1 <- get_pseudobulk(tmp$counts, clusters=tmp$celltypes, donors=tmp$donors)
  dim(test1) # should be 20 x 6
  test2 <- get_pseudobulk(tmp$counts, clusters=tmp$celltypes, donors=tmp$donors, trim=0)
  dim(test2) # should be 20 x 9
  sample <- sapply(strsplit(colnames(test2), "_"), function(x){x[[2]]})
  cell_type <- sapply(strsplit(colnames(test2), "_"), function(x){x[[1]]})
} 
\keyword{pseudobulk}
