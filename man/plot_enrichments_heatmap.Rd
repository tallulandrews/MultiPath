% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Visualization.R
\name{plot_enrichments_heatmap}
\alias{plot_enrichments_heatmap}
\title{Plot Enrichment Heatmap}
\usage{
plot_enrichments_heatmap(
  list_of_rich,
  pathways = NULL,
  pathways.names = NULL,
  ntop = 5,
  colors = (grDevices::colorRampPalette(rev(c("red", "orange", "yellow", "black", "navy",
    "purple", "magenta"))))(100),
  stars = TRUE,
  stars.col = "white",
  remove.prefix = TRUE,
  prefix.delim = "_",
  log.scale = FALSE,
  bounds = NULL,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  plot.result = TRUE,
  both.dir = FALSE,
  col.anno = NULL,
  row.anno = NULL
)
}
\arguments{
\item{list_of_rich}{a list of pathway enrichments for a single DE test obtained from do_ora or do_fgsea. All enrichments should be calculated using either do_ora or do_fgsea.}

\item{pathways}{a vector of pathway names of pathways to be plot.}

\item{pathways.names}{a vector of alias names for each pathway to be plot - used to shorted long pathway names when using predetermined pathways, overrides remove.prefix.}

\item{ntop}{number of pathways to select from each enrichement results.}

\item{colors}{a vector of colors for the heatmap}

\item{stars}{a boolean for whether to include stars to indicate significance.}

\item{stars.col}{the text colour of the stars}

\item{remove.prefix}{whether to remove the prefix from each pathway name.}

\item{prefix.delim}{the symbol that defines the end of the prefix}

\item{log.scale}{whether to log (base2) scale the enrichment scores}

\item{bounds}{boundarys for the colour scheme, if NULL then uses a centered boundaries for the maximum range of values.}

\item{cluster_rows}{whether to cluster the rows of the heatmap or not}

\item{cluster_cols}{whether to cluster the columns of the heatmap or not}

\item{plot.result}{whether to make the plot or just calculate the summarized results.}

\item{both.dir}{whether to select the ntop pathways both up & down or just the ntop based on pvalue regardless of direction}

\item{col.anno}{dataframe of annotations for each comparison. # Not implemented yet}

\item{row.anno}{dataframe of annotations for each pathway.  # Not implemented yet}
}
\value{
invisibly the data matrix and p value matrix used to make the heatmap (as a list)
}
\description{
Summarizes a set of pathway enrichments as a heatmap.
}
\details{
Selects the top pathways from each set of enrichment results. Then aggregates the significance and enrichment scores for each of the selected pathways across all enrichment results.
}
\examples{
list_of_rich <- lapply(rpois(5, lambda=10), generate_synthetic_enrichments, ngenes=100)
names(list_of_rich) <- paste("celltype", 1:length(list_of_rich), sep="")
plot_enrichments_heatmap(list_of_rich)
}
