#' Cell-type specific DE
#'
#' @description
#' Calculated differential expression for each cell-type given an appropriate design matrix
#'
#' @details
#' Uses edgeR (QLR test) to calculate differential expression for each coefficient in the design matrix for each cell-type. Results are filtered based on the provided FDR threshold.
#' 
#' @param pseudobulks a pseudobulk matrix created from get_pseudobulks
#' @param design_matrix a design matrix created with model.matrix
#' @param fdr False discovery rate threshold used to filter DE genes
#' @return a list of dataframes containing DE results for each cell-type
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' conditions <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[2]]})
#' names(conditions) <- colnames(example_celltype_pseudobulks)
#' design <- model.matrix(~conditions)
#' de <- compute_cell_type_specific_DE(example_celltype_pseudobulks, design)
#' @export
compute_cell_type_specific_DE <- function(pseudobulks, design_matrix, fdr=0.05) {
	rownames(design_matrix) <- colnames(pseudobulks)
	sample <- sapply(strsplit(colnames(pseudobulk), "_"), function(x){x[[2]]}) # correct order
	cell_type <- sapply(strsplit(colnames(pseudobulk), "_"), function(x){x[[1]]}) # correct order
	all_outs <- list()
	for (type in cell_type) {
		dat <- pseudobulks[,cell_type==type]
		these_samples <- sample[cell_type==type]
		# Create the design matrix for this cell-type
		design <- design_matrix[rownames(design_matrix) %in% these_samples,]
		if (matrixcalc::is.singular.matrix(design)) {
			warning(paste("Warning: DE for", type, "could not be computed because the predictors are dependent."))
			next;
		}
		all_outs[[type]] <- one_cell_type_DE(dat, design, fdr=fdr)
	}
	return(all_outs)
}

#' Cell-type specific DE - parallel
#'
#' @description
#' Calculated differential expression for each cell-type given an appropriate design matrix
#'
#' @details
#' Uses edgeR (QLR test) to calculate differential expression for each coefficient in the design matrix for each cell-type. Results are filtered based on the provided FDR threshold. Performs differential expression for each cell-type in parallel to speed up computations.
#' 
#' @param pseudobulks a pseudobulk matrix created from get_pseudobulks
#' @param design_matrix a design matrix created with model.matrix
#' @param fdr False discovery rate threshold used to filter DE genes
#' @param n.cores Number of CPUs to use for parallel computing
#' @return a list of dataframes containing DE results for each cell-type
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' conditions <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[2]]})
#' names(conditions) <- colnames(example_celltype_pseudobulks)
#' design <- model.matrix(~conditions)
#' de <- compute_cell_type_specific_DE_parallel(example_celltype_pseudobulks, design)
#' @export
compute_cell_type_specific_DE_parallel <- function(pseudobulks, design_matrix, fdr=0.05) {
	require(foreach)
	sample <- sapply(strsplit(colnames(pseudobulk_expr_mat_all), "_"), function(x){x[[1]]})
	cell_type <- sapply(strsplit(colnames(pseudobulk_expr_mat_all), "_"), function(x){x[[2]]})
	all_outs <- foreach( type=cell_type) %do% {
		dat <- pseudobulks[,cell_type==type]
		these_samples <- sample[cell_type==type]
		# Create the design matrix for this cell-type
		design <- design_matrix[rownames(design_matrix) %in% these_samples,]
		if (matrixcalc::is.singular.matrix(design)) {
			warning(paste("Warning: DE for", type, "could not be computed because the predictors are dependent."))
		} else {
			list(type=type, de=one_cell_type_DE(dat, design, fdr=fdr))
		}
	}
	type_names <- sapply(all_outs, function(x){x$type})
	all_outs <- lapply(all_outs, function(x){x$de})
	names(all_outs) <- type_names;
	return(all_outs)
}

#' Cell-type specific DE - one cell-type
#'
#' @description
#' Calculated differential expression for a cell-type given an appropriate design matrix
#'
#' @details
#' Called by compute_cell_type_specific_DE to calculate DE of each cell-type. This function filters genes for those with total counts > 5 and detection in at least 3 samples. Then runs DE using edgeR on the provided data & design matrix. Normalization is performed using TMM, and DE is tested using quasi-likelihood (QLF test). Each coefficient is tested and results are combined into a single dataframe. 
#' 
#' @param dat a pseudobulk matrix subset to a single cell-type
#' @param design a design matrix created with model.matrix for columns of dat
#' @param fdr False discovery rate threshold used to filter DE genes
#' @return a dataframe of DE results for each coefficient in the design matrix.
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' conditions <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[2]]})
#' cell_type <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[1]]})
#' names(conditions) <- colnames(example_celltype_pseudobulks)
#' conditions <- conditions[cell_type == cell_type[1]]
#' dat <- example_celltype_pseudobulks[,cell_type == cell_type[1]]
#' design <- model.matrix(~conditions)
#' de <- one_cell_type_DE(dat, design)
#' @export
one_cell_type_DE <- function(dat, design, fdr=0.05) {
	# Filter out lowly expressed genes
	dat <- dat[Matrix::rowSums(dat) > 5 & Matrix::rowSums(dat > 0) >= 3,]
	# Perform DE using edgeR
	dge <- edgeR::DGEList(counts=dat)
	dge <- calcNormFactors(dge, method="TMM", refColumn=which(colSums(dat) == max(colSums(dat))))
	dge <- estimateDisp(dge, design)
	fit <- glmQLFit(dge, design)
	all_contrast_outs <- c()
	for (i in 2:ncol(design)) {
		contrast_vec <- rep(0, ncol(design))
		contrast_vec[i] <- 1
		de <- topTags(dlmQLFTest(fit, contrast=contrast_vec, n=nrow(dat), p.value=fdr))
		de <- de$table
		de$Coefficient <- colnames(design)[i]
		de$Gene <- rownames(de)
		all_contrast_outs <- rbind(all_contrast_outs, de)
	}
	return(all_contrast_outs)
}
