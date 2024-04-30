#' rowMeans wrapper
#'
#' @description
#' A wrapper for Matrix::rowMeans that automatically checks if the input is a single column if so it returns just that column.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return the row means of the matrix or the original vector if supplied with a vector
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_rowMeans(example_data$counts)
#' my_rowMeans(c(1,2,3,4,5,6,7,8,10))
#' @export
my_rowMeans <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowMeans(x))
                }
        }
        return(x);
}

#' rowSums wrapper
#'
#' @description
#' A wrapper for Matrix::rowSums that automatically checks if the input is a single column if so it returns just that column.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return the row sums of the matrix or the original vector if supplied with a vector
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_rowSums(example_data$counts)
#' my_rowSums(c(1,2,3,4,5,6,7,8,10))
#' @export
my_rowSums <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowSums(x))
                }
        }
        return(x);
}

#' colMeans wrapper
#'
#' @description
#' A wrapper for Matrix::colMeans that automatically checks if the input is a single row if so it returns just that row.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return the column means of the matrix or the original vector if supplied with a vector
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_colMeans(example_data$counts)
#' my_colMeans(c(1,2,3,4,5,6,7,8,10))
#' @export
my_colMeans <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colMeans(x))
                }
        }
        return(x);
}

#' colSums wrapper
#'
#' @description
#' A wrapper for Matrix::colSums that automatically checks if the input is a single row if so it returns just that row.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return the column means of the matrix or the original vector if supplied with a vector
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_colSums(example_data$counts)
#' my_colSums(c(1,2,3,4,5,6,7,8,10))
#' @export
my_colSums <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colSums(x))
                }
        }
        return(x);
}

#' Sparse-Matrix compatible rowVars
#'
#' @description
#' Uses rowMeans to calculate variance of each row for sparse matrices.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return variances of each row of the matrix or 0s if supplied with a vector.
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_rowVars(example_data$counts)
#' my_rowVars(c(1,2,3,4,5,6,7,8,10))
#' @export
my_rowVars <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
			center <- Matrix::rowMeans(x)
			n <- ncol(x)
			vars <- n/(n-1)*(Matrix::rowMeans(x^2) - center^2)
                        return(vars)
                }
        }
        return(rep(0, length(x)));
}

#' Sparse-Matrix compatible colVars
#'
#' @description
#' Uses colMeans to calculate variance of each column for sparse matrices.
#' 
#' @param x a matrix or Matrix or a vector.
#' @return variances of each column of the matrix or 0s if supplied with a vector.
#' @examples
#' example_data <- generate_test_cellcounts()
#' my_colVars(example_data$counts)
#' my_colVars(c(1,2,3,4,5,6,7,8,10))
#' @export
my_colVars <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
			center <- Matrix::colMeans(x)
			n <- nrow(x)
			vars <- n/(n-1)*(Matrix::colMeans(x^2) - center^2)
                        return(vars)
                }
        }
        return(rep(0, length(x)));
}

#' Group Row Means
#'
#' @description
#' Optimized code for calculating the row means or row sums of a matrix after grouping the columns by a vector of groups.
#' 
#' @details
#' Uses my_rowSums, my_rowMeans, and my_rowVars to calculate a matrix of row means, row sums, or row variances for each group devinfed in group_labs. Note normalized expression should be averaged while raw counts should be summed.
#' @param MAT a matrix or sparse matrix of values to be aggregated
#' @param group_labs a vector of group labels equal in length to the number of columsn of MAT.
#' @param type whether to calculate the row means or row sums
#' @return a matrix of the row means for each group.
#' @examples
#' example_data <- generate_test_cellcounts()
#' total_celltype_counts <- group_rowmeans(example_data$counts,example_data$celltypes, type="sum")
#' avg_donor_counts <- group_rowmeans(example_data$counts,example_data$donors, type="mean")
#' #Test cases
#' group_rowmeans(example_data$counts, c("A", rep("B", ncol(example_data$counts)-1)), type="sum") # group with 1 entry
#' group_rowmeans(example_data$counts, c("A", rep("B", ncol(example_data$counts)-1)), type="mean") # group with 1 entry
#' group_rowmeans(example_data$counts, c("A", rep("B", ncol(example_data$counts)-1)), type="var") # group with 1 entry
#' group_rowmeans(example_data$counts[1,], example_data$donor, type="sum") # data with 1 row
#' group_rowmeans(example_data$counts[1,], example_data$donor, type="mean") # data with 1 row
#' group_rowmeans(example_data$counts[1,], example_data$donor, type="var") # data with 1 row
#' @export
group_rowmeans <- function(MAT, group_labs, type=c("mean","sum", "var")) {
        d <- split(seq(ncol(MAT)), group_labs);
	if (type[1] == "mean") {
		if(nrow(MAT) > 1) {
	        	mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
		} else {
			mus <- sapply(d, function(group) mean(MAT[,group])) # only one row
		}
	} 
	if (type[1] == "sum") {
		if (nrow(MAT) > 1) {
	        	mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
		} else {
			mus <- sapply(d, function(group) sum(MAT[,group])) # only one row
		}
	} 
	if (type[1] == "var") {
		if (nrow(MAT) > 1) {
	        	mus <- sapply(d, function(group) my_rowVars(MAT[,group]))
		} else {
			mus <- sapply(d, function(group) stats::var(MAT[,group])) # only one row
		}
	}
        return(mus);
}

#' Group Column Means
#'
#' @description
#' Optimized code for calculating the column means or column sums of a matrix after grouping the columns by a vector of groups.
#' 
#' @details
#' Uses my_colSums, my_colMeans, and my_rowVars to calculate a matrix of column means, column sums, or column variances for each group devinfed in group_labs. Note normalized expression should be averaged while raw counts should be summed.
#' @param MAT a matrix or sparse matrix of values to be aggregated
#' @param group_labs a vector of group labels equal in length to the number of columsn of MAT.
#' @param type whether to calculate the column means, column sums, or column variances
#' @return a matrix of the column means for each group.
#' @examples
#' example_data <- generate_test_cellcounts()
#' total_celltype_counts <- group_colmeans(t(example_data$counts),example_data$celltypes, type="sum")
#' avg_donor_counts <- group_colmeans(t(example_data$counts),example_data$donors, type="mean")
#' #Test cases
#' group_colmeans(example_data$counts, c("A", rep("B", nrow(example_data$counts)-1)), type="sum") # group with 1 entry
#' group_colmeans(example_data$counts, c("A", rep("B", nrow(example_data$counts)-1)), type="mean") # group with 1 entry
#' group_colmeans(example_data$counts, c("A", rep("B", nrow(example_data$counts)-1)), type="var") # group with 1 entry
#' group_colmeans(t(example_data$counts)[,1], example_data$donor, type="sum") # data with 1 column
#' group_colmeans(t(example_data$counts)[,1], example_data$donor, type="mean") # data with 1 column
#' group_colmeans(t(example_data$counts)[,1], example_data$donor, type="var") # data with 1 column
#' @export
group_colmeans <- function(MAT, group_labs, type=c("mean", "sum", "var")) {
        d <- split(seq(nrow(MAT)), group_labs);
	if (type[1] == "mean") {
		if(ncol(MAT) > 1) {
        		mus <- sapply(d, function(group) my_colMeans(MAT[group,]))
		} else {
			mus <- sapply(d, function(group) mean(MAT[group,])) # only one col
		}
	}
	if (type[1] == "sum") {
		if(ncol(MAT) > 1) {
        		mus <- sapply(d, function(group) my_colSums(MAT[group,]))
		} else {
			mus <- sapply(d, function(group) sum(MAT[group,])) # only one col
		}
	}
	if (type[1] == "var") {
		if(ncol(MAT) > 1) {
        		mus <- sapply(d, function(group) my_colVars(MAT[group,]))
		} else {
			mus <- sapply(d, function(group) stats::var(MAT[group,])) # only one col
		}
	}
        return(mus);
}


#' Trim Pseudobulk Samples
#'
#' @description
#' Removes sample x cell-type combinations where there are too few cells to obtain reliable pseudobulk estimates.
#' 
#' @details
#' Any combination of cell-type x individual (or sample) that has fewer than `nmin` cells are identified for removal. This is called by get_pseudobulks.
#' @param clusters a vector of cell-type names or cluster IDs.
#' @param individual a vector of individual or sample IDs.
#' @param nmin the minimum number of cells required for a reliable pseudobulk estimate.
#' @return a vector of cluster-individual pairs that should be excluded
#' @examples
#' example_data <- generate_test_cellcounts()
#' trim_for_pseudobulks(example_data$celltypes, example_data$donors) # all of celltype1 should be identified for removal.
#' @export
trim_for_pseudobulks <- function(clusters, individual, nmin=10) {
        tmp <- table(clusters, individual)
        toofew <- which(tmp < nmin, arr.ind=TRUE)
        to_exclude <- paste(rownames(tmp)[toofew[,1]], colnames(tmp)[toofew[,2]])
        exclude <- paste(clusters, individual) %in% to_exclude
        return(exclude)
}

#' Generate Pseudobulks
#'
#' @description
#' Calculates pseudobulk expression for each celltype x sample pair.
#' 
#' @details
#' Calculates pseudobulk expression for performing differential expression (DE). If using a negative-binomial model based DE method such as edgeR or DESeq2 the "sum" of the raw umi counts should be used. If using a Gaussian (e.g. MAST) or non-parametric DE method then the "mean" of normalized expression should be calculated.
#' To avoid noise from samples where only a small number of cells of a particular cell-type are found, it automatically filters out cases where there are fewer than `trim` cells of a particular cell-type in a particular sample.
#' Note: to avoid issues retreiving sample & cell-type IDs from the output column names, all underscores in the original clusters and donor ids are replaced with a dash("-").
#' 
#' @param mat a matrix of umi counts or normalized expression for each cell, genes=rows, cells=columns, supports sparse matrices
#' @param clusters a vector of cluster or cell-type labels for each cell.
#' @param individual a vector of patient or sample labels for each cell.
#' @param method whether to add up or average the expression in each cell-type x sample group of cell.
#' @param trim cell-type x sample groups with fewer than this many cells will not be included in the pseudobulk matrix
#' @param refactor whether to factor the clusters & donors vectors - recommended if trimming.
#' @return  A matrix of pseudobulk expression with column names in the format: [cluster]_[donor].
#' @examples
#'   example_data <- generate_test_cellcounts()
#'   test1 <- get_pseudobulk(example_data$counts, clusters=example_data$celltypes, donors=example_data$donors)
#'   dim(test1) # should be 20 x 6
#'   test2 <- get_pseudobulk(example_data$counts, clusters=example_data$celltypes, donors=example_data$donors, trim=0)
#'   dim(test2) # should be 20 x 9
#'   sample <- sapply(strsplit(colnames(test2), "_"), function(x){x[[2]]})
#'   cell_type <- sapply(strsplit(colnames(test2), "_"), function(x){x[[1]]})
#' @export
get_pseudobulk <- function(mat, clusters, individual, method=c("sum", "mean"), trim=10, refactor=TRUE) {
	#avoid naming issues - we use underscores to sepearate cluster vs donor in pseudobulk column names
	clusters <- sub("_", "-", clusters)
	individual <- sub("_", "-", individual)

	# Remove cluster-donor pairs where there are too few cells to get a reliable expression profile
	if (trim > 0) {
		to.exclude <- trim_for_pseudobulks(clusters, individual, nmin=trim)
		mat <- mat[,!to.exclude]
		clusters <- clusters[!to.exclude]
		individual <- individual[!to.exclude]
	}
	if (refactor) {
		clusters <- factor(clusters)
		individual <- factor (individual)
	}

	# Subset "mat" to each pair of donors & clusters and add up the umi counts for all the cells
        c <- split(seq(ncol(mat)), clusters);
        # expression per donor in this cluster
        clust_expr <- lapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], individual[clust], type=method[1]);
                if(is.null(dim(d_expr))) {
                        l <- sapply(d_expr, length)
                        keep <- which(l == nrow(mat))
                        d_expr <- matrix(unlist(d_expr[keep]), ncol=length(keep), byrow=FALSE);
                        rownames(d_expr) <- rownames(mat);
                        colnames(d_expr) <- paste(clusters[clust[1]], levels(individual)[keep], sep="_")
                } else {
                        colnames(d_expr) <- paste(clusters[clust[1]], colnames(d_expr), sep="_")
                }
                return(d_expr);
        })
        out <- clust_expr[[1]];
	c_names <- colnames(out)
        for (i in 2:length(clust_expr)) {
                c_names <- c(c_names, colnames(clust_expr[[i]]))
                out <- cbind(out, clust_expr[[i]]);
                if (is.null(dim(out))){
                        out <- matrix(out, ncol=1)
                        rownames(out) <- rownames(mat)
                }
                colnames(out) <- c_names
        }
        return(out)
}

#' Remove Duplicate Rows
#'
#' @description
#' Removes rows from a matrix where there are duplicate row names, e.g. when changing gene IDs between databases or species.
#' 
#' @details
#' Combines or selects one row of a matrix to represent all rows with duplicated rownames. There are three approaches:
#' max : keep the row with the highest average across columns
#' sum : add up the rows with duplicated rownames for each column.
#' mean: take the average across the duplicated rows for each column 
#' @param mat a matrix or sparse matrix of numeric data, such as a gene expression matrix.
#' @param rownames_dups a vector of rownames for `mat` that contains duplicates.
#' @param method which approach to use to aggregate data across duplicated rows (see:Details).
#' @return a matrix with the provided rownames but without duplicates.
#' @examples
#' example_data <- generate_test_cellcounts()$counts
#' gene_ids <- sample(paste("gene", 1:10, sep=""), 20, replace=TRUE)
#' deduped <- remove_duplicate_rows(example_data, gene_ids, method="max")
#' total <- remove_duplicate_rows(example_data, gene_ids, method="sum")
#' avg <- remove_duplicate_rows(example_data, gene_ids, method="mean")
#' dim(deduped) == dim(total) # TRUE
#' dim(deduped) == dim(avg) # TRUE
#' dim(deduped)[1] == length(unique(gene_ids)) # TRUE
#' @export
remove_duplicate_rows <- function(mat, rownames_dups, method=c("max", "sum", "mean")) {
        rownames_dups <- as.character(rownames_dups)
        for (g in unique(rownames_dups[duplicated(rownames_dups)])) {
                print(paste(g, "is duplicated", sum(rownames_dups == g), "times"))
                b_means <- rowMeans(mat)
                to.remove <- which(rownames_dups == g)
                if (method[1] == "max") {
                        top <- max(b_means[to.remove])
                        to.remove <- rownames_dups == g & b_means < top
                } else {
                        if (method[1] == "sum") {
                                new_row <- colSums(mat[to.remove,])
                        } else if (method[1] == "mean") {
                                new_row <- colMeans(mat[to.remove,])
                        }
                        mat <- rbind(mat, new_row)
                        rownames_dups <- c(rownames_dups, g)
                }
                mat <- mat[!to.remove,]
                rownames_dups <- rownames_dups[!to.remove]
        }
        final_is.dup <- duplicated(rownames_dups);
        mat <- mat[!final_is.dup,]
        rownames_dups <- rownames_dups[!final_is.dup]
        rownames(mat) <- rownames_dups
        return(mat)
}
