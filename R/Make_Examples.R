#' Generate synthetic pseudobulks (one-celltype)
#'
#' @description
#' Generates a synthetic pseudobulk matrix for testing purposes.
#' 
#' @details
#' Generates synthetic pseudobulk data for one cell-type. Sample-specific library sizes are drawn from an exponential distribution and rescaled to the median value. Gene-specific means are drawn from an exponential distribution with rate=1/100, counts are drawn from a negative binomial distribution.
#' @param ngenes number of genes to generate data for.
#' @param de_prop proportion of genes to be differentially expressed between biological conditions.
#' @param nconditions number of different biological conditions for which to generate synthetic data.
#' @param nsamples_per_condition number of samples for each biological condition.
#' @param dispersion_factor dispersion parameter (1/size) of the negative bionomial distribution.
#' @return A matrix of pseudobulk expression with columns named as: [celltype]_[condition]-[sample]
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks_one_cell_type()
#' dim(example_celltype_pseudobulks)
#' @export
generate_synthetic_pseudobulks_one_cell_type <- function(ngenes=100, de_prop=0.5, nconditions=2, nsamples_per_condition=3, dispersion_factor=0.2){
	if (nconditions < 2) {stop("Must have at least 2 biological conditions.")}
	nsamples = nsamples_per_condition*nconditions

	lib_size <- stats::rexp(nsamples, rate=1)
	lib_size <- lib_size/stats::median(lib_size)
	g_means <- stats::rexp(ngenes, rate=1/100)
	g_mean_mat <- matrix(rep(g_means, nsamples_per_condition*nconditions), ncol=nsamples_per_condition*nconditions, nrow=ngenes)

	de <- 1:round(ngenes*de_prop);

	for (i in 1:nconditions) {
		these_samples = ((i-1)*nsamples_per_condition+1):((i-1)*nsamples_per_condition+nsamples_per_condition)
		these_g_means <- stats::rexp(max(de), rate=1/100)
		g_mean_mat[de, these_samples] <- matrix(rep(these_g_means, nsamples_per_condition), ncol=nsamples_per_condition, nrow=max(de))
	}
	g_mean_mat <- t(round(t(g_mean_mat)*lib_size))

	expr_counts <- t(apply(g_mean_mat, 1, function(x) {
		sapply(x, function(y){
			stats::rnbinom(1, mu=y, size=1/dispersion_factor)
			})
		}))
	rownames(expr_counts) <-  paste("Gene", 1:ngenes, sep="")
	colnames(expr_counts) <-  paste("C", rep(1:nconditions, each=nsamples_per_condition),"-I",1:(nconditions*nsamples_per_condition), sep="")
	return(expr_counts)
}

#' Generate synthetic pseudobulks
#'
#' @description
#' Generates a synthetic pseudobulk matrix for testing purposes.
#' 
#' @details
#' Generates synthetic pseudobulk data for a whole matrix. Sample-specific library sizes are drawn from an exponential distribution and rescaled to the median value. Gene-specific means are drawn from an exponential distribution with rate=1/100, counts are drawn from a negative binomial distribution.
#' @param ngenes number of genes to generate data for.
#' @param nconditions number of different biological conditions for which to generate synthetic data.
#' @param nsamples_per_condition number of samples for each biological condition.
#' @param ncell_types number of cell-types.
#' @param dispersion_factor dispersion parameter (1/size) of the negative bionomial distribution.
#' @return A matrix of pseudobulk expression with columns named as: [celltype]_[condition]-[sample]
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' dim(example_celltype_pseudobulks)
#' @export

generate_synthetic_pseudobulks <- function(ngenes=100, nconditions=2, nsamples_per_condition=3, ncell_types=3, dispersion_factor=0.2) {
	all_pseudobulks <- c();
	for (type in 1:ncell_types) {
		this_type <- paste("celltype", type, sep="-")
		this_type_counts <- generate_synthetic_pseudobulks_one_cell_type(ngenes=100, nconditions=2, nsamples_per_condition=3, dispersion_factor=0.2)
		colnames(this_type_counts) <- paste(paste("celltype", type, sep=""), colnames(this_type_counts), sep="_")
		all_pseudobulks <- cbind(all_pseudobulks, this_type_counts)
	}
	return(all_pseudobulks)
}

## Excessively Complicated
#generate_synthetic_cellcounts <- function(ngenes=100, nclusters=10, avg_ncells_per_donor=100, ndonor=6){
#	g_means <- rexp(ngenes, rate=1/100)
#	cells_per_cluster <- rpois(ndonor*nclusters, lambda=avg_ncells_per_donor/nclusters)
#	total_cells <- sum(cells_per_cluster)
#	
#	cells_per_donor <- sapply(split(1:length(cells_per_cluster),rep(1:nclusters, each=ndonor)), function(x){sum(cells_per_cluster[x])})
#	donor_lab <- rep(paste("donor", 1:ndonor, sep=""), each=cells_per_donor)
#	
#	mat <- matrix(rnorm(1000), ncol=100)

#' Generate Cell Counts for testing
#'
#' @description
#' Generates a synthetic gene x cell umi count matrix for testing methods.
#' 
#' @details
#' Generates a matrix of uniform Poisson distributed counts, as well as vectors of cell-type
#' and donor of origin labels. These data include one cell-type represented by 1 cell per donor, 
#' one with 15 cells/donor and one with 30 cells/donor to check for errors in methods.
#' @return a list with items: counts = umi count matrix, donors = vector of sample of origin labels, celltype = vector of celltype labels
#' @examples
#' example_data <- generate_test_cellcounts()
#' dim(example_data$counts)
#' table(example_data$donors)
#' table(example_data$celltypes)
#' @export

generate_test_cellcounts <- function() {
	ndonor=3
	ngene=20;
	celltypes <- c(rep("type1", ndonor), rep("type2", ndonor*15), rep("type3", ndonor*30))
	donors <- rep(c("donor1", "donor2", "donor3"), times=(1+15+30))
	mat <- matrix(stats::rpois(ngene*length(celltypes), lambda=10), ncol=length(celltypes))
	colnames(mat) <- paste("cell", 1:ncol(mat), sep="")
	rownames(mat) <- paste("gene", 1:nrow(mat), sep="")
	return(list(counts=mat, donors=donors, celltypes=celltypes))
}

#' Generate synthetic enrichments
#'
#' @description
#' Generates a synthetic enrichments for testing other functions
#' 
#' @details
#' Generates a set of random output formatted to look like the output from do_ora or do_fgsea.
#' @param npathways number of pathway enrichments to generate
#' @param ngenes number of total genes to use to populate enrichments.
#' @return A list containing a dataframe of the same format as the output of do_ora or do_fgsea and a list of contributing genes.
#' @examples
#' example_enrichments <- generate_synthetic_enrichments(10)
#' dim(example_enrichments$results)
#' length(example_enrichments$contrib)
#' @export

generate_synthetic_enrichments <- function(npathways=10, ngenes=200){
	all_genes <- paste("Gene", 1:ngenes, sep="")
	contrib <- list() # genes contributing to each pathway
	for (i in 1:npathways) {
		contrib[[i]] <- unique(sample(all_genes, size=stats::rpois(1, lambda=12), replace=TRUE))
	}
	intersection <- sapply(contrib, length)
	res <- data.frame(pathway=paste("pathway",1:npathways,sep=""),
			intersection=intersection,
			log2fe=stats::rnorm(npathways, sd=2),
			FDR=exp(-abs(stats::rnorm(npathways,mean=5, sd=10))))
	rownames(res) <- res$pathway
	names(contrib) <- res[,1]
	return(list(results=res, contrib=contrib))
}

## Excessively Complicated
#generate_synthetic_cellcounts <- function(ngenes=100, nclusters=10, avg_ncells_per_donor=100, ndonor=6){
#	g_means <- rexp(ngenes, rate=1/100)
