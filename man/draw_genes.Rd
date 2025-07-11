% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hierarchical_Pathways.R
\name{draw_genes}
\alias{draw_genes}
\title{Draw Genes from Pathways}
\usage{
draw_genes(gt_pathways, N = 50, fpr = 0.05, genes, gene_anno, seed = NULL)
}
\arguments{
\item{gt_pathways}{a vector of pathways to draw genes from.}

\item{N}{the number of genes to draw from each pathway.}

\item{fpr}{the proportion of the returned genes that should be false-positives.}

\item{genes}{a vector of gene IDs that constitute the background - i.e. all detectable genes.}

\item{gene_anno}{a table of gene-pathway annotations, must have a column named "GO" containing the pathway term, and a column named "GENE" with the gene ids.}

\item{seed}{the random seed.}
}
\value{
a vector of selected genes.
}
\description{
Randomly draws true positive and false positive genes from a set of provided pathways.
}
\details{
Given a set of pathways and a set of tissue-specific genes, this function randomly selects a set number of genes from each pathway that are also tissue-specific. Subsequently randomly selects genes uniformly from the tissue-specific genes as false-positive genes. Warning: no error checking / input validation is preformed for this function.
}
\examples{
gene2go <- load_GOBP_info()
tissue_specific <- get_tissue_specific_genes(1, seed=301)
genes <- draw_genes(c("GO:1990573","GO:0051898"), N=5, fpr=0.5, tissue_specific[[1]], gene2go, seed=101)
genes <- draw_genes(c("GO:1990573","GO:0051898"), N=5, fpr=0, tissue_specific[[1]], gene2go, seed=101)
}
