my_rowMeans <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowMeans(x))
                }
        }
        return(x);
}
my_rowSums <- function(x) {
        if (!is.null(ncol(x))) {
                if (ncol(x) > 1) {
                        return(Matrix::rowSums(x))
                }
        }
        return(x);
}
my_colMeans <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colMeans(x))
                }
        }
        return(x);
}
my_colSums <- function(x) {
        if (!is.null(nrow(x))) {
                if (nrow(x) > 1) {
                        return(Matrix::colSums(x))
                }
        }
        return(x);
}

group_rowmeans <- function(MAT, group_labs, type=c("mean","sum")) {
        d <- split(seq(ncol(MAT)), group_labs);
	if (type[1] == "mean") {
	        mus <- sapply(d, function(group) my_rowMeans(MAT[,group]))
	} else {
	        mus <- sapply(d, function(group) my_rowSums(MAT[,group]))
	} 
        return(mus);
}
group_colmeans <- function(MAT, group_labs, type=c("mean", "sum")) {
        d <- split(seq(nrow(MAT)), group_labs);
	if (type[1] == "mean") {
        	mus <- sapply(d, function(group) my_colMeans(MAT[group,]))
	} else {
        	mus <- sapply(d, function(group) my_colSums(MAT[group,]))
	}
        return(mus);
}


trim_for_pseudobulks <- function(clusters, individual, nmin=10) {
        tmp <- table(clusters, individual)
        toofew <- which(tmp < nmin, arr.ind=TRUE)
        to_exclude <- paste(rownames(tmp)[toofew[,1]], colnames(tmp)[toofew[,2]])
        exclude <- paste(clusters, individual) %in% to_exclude
        return(exclude)
}

# Table of total expression of cells from each donor in each cluster
#  - for edgeR
# Calculate pseudobulks for each cluster in each individual for scRNAseq experiments with multiple replicates.
# For edgeR or DESeq2 we recommend using method="sum" and the raw counts as the matrix.
# if using a standard glm or MAST we recommend using method="mean" on the lognormalized data
# trim provides the option to remove groups where there are insufficient samples to get a reliable estimate of the mean (i.e. fewer than trim cells). Set trim <0 to turn off trimming.
# refactor determines whether clusters & donors are refactored. Refactoring will change the order of the groups in the output but will avoid errors where some clusters contain 0 cells. It is highly recommended to refactor if using trim.
get_pseudobulk <- function(mat, clusters, donors, method=c("sum", "mean"), trim=10, refactor=TRUE) {
	#avoid naming issues - we use underscores to sepearate cluster vs donor in pseudobulk column names
	clusters <- sub("_", "-", clusters)
	donor <- sub("_", "-", donors)

	# Remove cluster-donor pairs where there are too few cells to get a reliable expression profile
	if (trim > 0) {
		to.exclude <- trim_for_pseudobulks(clusters, donors, nmin=trim)
		mat <- mat[,!to.exclude]
		clusters <- clusters[!to.exclude]
		donors <- donors[!to.exclude]
	}
	if (refactor) {
		clusters <- factor(clusters)
		donors <- factor (donors)
	}

	# Subset "mat" to each pair of donors & clusters and add up the umi counts for all the cells
        c <- split(seq(ncol(mat)), clusters);
        donor_freqs <- table(donors)/length(donors)
        # expression per donor in this cluster
        clust_expr <- sapply(c, function(clust) {
                d_expr <- group_rowmeans(mat[,clust], donors[clust], type=method[1]);
                if(is.null(dim(d_expr))) {
                        l <- sapply(d_expr, length)
                        keep <- which(l == nrow(mat))
                        d_expr <- matrix(d_expr[[keep]], ncol=length(keep), byrow=FALSE);
                        rownames(d_expr) <- rownames(mat);
                        colnames(d_expr) <- paste(clusters[clust[1]], levels(donors)[keep], sep="_")
                } else {
                        colnames(d_expr) <- paste(clusters[clust[1]], colnames(d_expr), sep="_")
                }
                return(d_expr);
        })
        out <- clust_expr[[1]];
        for (i in 2:length(clust_expr)) {
                c_names <- c(colnames(out), colnames(clust_expr[[i]]))
                out <- cbind(out, clust_expr[[i]]);
                if (is.null(dim(out))){
                        out <- matrix(out, ncol=1)
                        rownames(out) <- rownames(mat)
                }
                colnames(out) <- c_names
        }
        return(out)
}

# remove duplicate rows, this can be handy if remapping gene IDs.
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
