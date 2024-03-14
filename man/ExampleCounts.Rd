\name{generate_test_cellcounts}
\alias{generate_test_cellcounts}
\title{Generate Exmple Counts}
\description{
	Generates a random count matrix for unit testing.
}
\usage{
	generate_test_cellcounts()
}
\arguments{
}
\details{ 
	Generates a random matrix (poisson distributed) containing 20 genes for 3 cell-types and 3 donors where each donor has 1 cell of type1, 15 cells of type2 and 30 cells of type3.
}
\value{
	A random count matrix.
}
\examples{
	generate_test_cellcounts()
} 
\keyword{unit-test}
