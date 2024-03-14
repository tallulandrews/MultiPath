\name{generate_synthetic_pseudobulks}
\alias{generate_synthetic_pseudobulks}
\alias{generate_synthetic_pseudobulks_one_cell_type}
\title{Generate Example Pseudobulks}
\description{
	Generates a synthetic pseudobulk matrix for testing purposes.
}
\usage{
	generate_synthetic_pseudobulks(ngenes=100, nconditions=2, nsamples_per_condition=3, ncell_types=3, dispersion_factor=0.2)
	generate_synthetic_pseudobulks_one_cell_type(ngenes=100, de_prop=0.5, nconditions=2, nsamples_per_condition=3, dispersion_factor=0.2)
}
\arguments{
  \item{ngenes}{number of genes to generate data for.}
  \item{de_prop}{proportion of genes to be differentially expressed between biological conditions}
  \item{nconditions}{number of different biological conditions to generate synthetic datafor}

  \item{nsamples_per_condition}{number of samples for each biological condition.}
  \item{ncell_types}{number of cell-types.}
  \item{dispersion_factor}{dispersion parameter (1/size) of the negative bionomial distribution.}
}
\details{ 
	Generates synthetic pseudobulk data for a whole matrix or one cell-type at a time. Sample-specific library sizes are drawn from an exponential distribution and rescaled to the median value. Gene-specific means are drawn from an exponential distribution with rate=1/100, counts are drawn from a negative binomial distribution.
}
\value{
	A matrix of pseudobulk expression with columns named as: [celltype]_[condition]-[sample]
}
\examples{
  example_celltype_pseudobulks <- generate_synthetic_pseudobulks_one_cell_type()
  example_data <- generate_synthetic_pseudobulks()
} 
\keyword{pseudobulk}
