generate_synthetic_pseudobulks_one_cell_type <- function(ngenes=100, de_prop=0.5, nconditions=2, nsamples_per_condition=3, dispersion_factor=0.2){
	if (nsamples < nconditions) {stop("Must have more samples than conditions.")}
	if (nconditions < 2) {stop("Must have at least 2 biological conditions.")}

	lib_size <- rexp(nsamples, rate=1)
	lib_size <- lib_size/median(lib_size)
	g_means <- rexp(ngenes, rate=1/100)
	g_mean_mat <- matrix(rep(g_means, nsamples_per_condition*nconditions), ncol=nsamples_per_condition*nconditions, nrow=ngenes)

	de <- 1:round(ngenes*de_prop);

	for (i in 1:nconditions) {
		these_samples = ((i-1)*samples_per_condition+1):((i-1)*samples_per_condition+samples_per_condition)
		these_g_means <- rexp(max(de), rate=1/100)
		g_mean_mat[de, these_samples] <- matrix(rep(these_g_means, samples_per_condition), ncol=samples_per_condition, nrow=max(de))
	}
	g_mean_mat <- t(round(t(g_mean_mat)*lib_size))

	expr_counts <- t(apply(g_mean_mat, 1, function(x) {
		sapply(x, function(y){
			rnbinom(1, mu=y, size=1/dispersion_factor)
			})
		}))
	rownames(expr_counts) <-  paste("Gene", 1:ngenes, sep="")
	colnames(expr_counts) <-  paste("C", rep(1:nconditions, each=nsamples_per_condition),"-I",1:(nconditions*nsamples_per_condition), sep="")
	return(expr_counts)
}

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

generate_test_cellcounts <- function() {
	ndonor=3
	ngene=20;
	celltypes <- c(rep("type1", ndonor), rep("type2", ndonor*15), rep("type3", ndonor*30))
	donors <- rep(c("donor1", "donor2", "donor3"), times=(1+15+30))
	mat <- matrix(rpois(ngene*length(celltypes), lambda=10), ncol=length(celltypes))
	colnames(mat) <- paste("cell", 1:ncol(mat), sep="")
	rownames(mat) <- paste("gene", 1:nrow(mat), sep="")
	return(list(counts=mat, donors=donors, celltypes=celltypes))
}
