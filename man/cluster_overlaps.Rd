% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Enrichments.R
\name{cluster_overlaps}
\alias{cluster_overlaps}
\title{Cluster Pathways}
\usage{
cluster_overlaps(overlaps, equivalent = 0.5, plot.result = TRUE)
}
\arguments{
\item{overlaps}{a matrix of pair-wise overlaps between contributing genes to each pathway, this is obtained from 'get_overlaps'.}

\item{equivalent}{a scalar value for the equivalence threshold.}

\item{plot.result}{a boolean determining whether to plot a heatmaps of the overlaps showing he pathway groups.}
}
\value{
A vector of group IDs for each pathway.
}
\description{
Clusters pathways based on overlapping contributing genes.
}
\details{
Uses single-linkage hierarchical clustering(`hclust`) to group pathway terms based on overlapping contributing genes.
}
\examples{
#paths <- convert_GSEAObj_to_list(get_pathways("test"))
paths <- test_paths;print("2")
rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)));
overlaps <- get_overlaps(rich);
pathway_groups <- cluster_overlaps(overlaps)
}
