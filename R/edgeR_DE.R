#' Cell-type specific DE
#'
#' @description
#' Calculated differential expression for each cell-type given an appropriate design matrix
#'
#' @details
#' Uses edgeR (QLR test) or DESeq2 to calculate differential expression for each coefficient in the design matrix for each cell-type. Results are filtered based on the provided FDR threshold. See: one_cell_type_DE_edgeR and one_cell_type_DE_deseq2 for details.
#' 
#' @param pseudobulks a pseudobulk matrix created from get_pseudobulks
#' @param design_matrix a design matrix created with model.matrix
#' @param fdr False discovery rate threshold used to filter DE genes
#' @param de_method Which pseudobulk DE method to use (deseq2 or edger).
#' @return a list of dataframes containing DE results for each cell-type
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' conditions <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[2]]})
#' names(conditions) <- colnames(example_celltype_pseudobulks)
#' design <- model.matrix(~conditions)
#' de <- compute_cell_type_specific_DE(example_celltype_pseudobulks, design, de_method="deseq2")
#' de <- compute_cell_type_specific_DE(example_celltype_pseudobulks, design, de_method="edger")
#' @export
compute_cell_type_specific_DE <- function(pseudobulks, design_matrix, fdr=0.05, de_method=c("deseq2", "edger")) {
	rownames(design_matrix) <- colnames(pseudobulks)
	sample <- sapply(strsplit(colnames(pseudobulks), "_"), function(x){x[[2]]}) # correct order
	cell_type <- sapply(strsplit(colnames(pseudobulks), "_"), function(x){x[[1]]}) # correct order
	all_outs <- list()
	if (is.null(cell_types)) {
		cell.types <-  unique(cell_type)
	}
	for (type in cell.types) {
		print(type)
		if (sum(cell_type==type) < 2) {
			warning(paste("Warning: DE for", type, "could not be computed because fewer than 2 samples contain this type."))
			next;
		}
		dat <- pseudobulks[,cell_type==type]
		# Create the design matrix for this cell-type
		design <- design_matrix[rownames(design_matrix) %in% colnames(dat),]
		if (Matrix::rankMatrix(design)[1] < ncol(design)) {
			warning(paste("Warning: DE for", type, "could not be computed because the predictors are dependent."))
			next;
		}
		if(de_method[1] == "deseq2") {
			all_outs[[type]] <- one_cell_type_DE_edgeR(dat, design, fdr=fdr)
		} else if (de_method[1] == "edger") {
			all_outs[[type]] <- one_cell_type_DE_deseq2(dat, design, fdr=fdr)
		} else {
			print("Error: unrecognized DE method")
		}
	}
	return(all_outs)
}

#' Cell-type specific DE - parallel
#'
#' @description
#' Calculated differential expression for each cell-type given an appropriate design matrix
#'
#' @details
#' Uses edgeR (QLR test) or DESeq2 to calculate differential expression for each coefficient in the design matrix for each cell-type. Results are filtered based on the provided FDR threshold. Performs differential expression for each cell-type in parallel to speed up computations. See: one_cell_type_DE_edgeR and one_cell_type_DE_deseq2 for details.
#' 
#' @param pseudobulks a pseudobulk matrix created from get_pseudobulks
#' @param design_matrix a design matrix created with model.matrix
#' @param fdr False discovery rate threshold used to filter DE genes
#' @param n.cores Number of CPUs to use for parallel computing
#' @param de_method Which pseudobulk DE method to use (deseq2 or edger).
#' @return a list of dataframes containing DE results for each cell-type
#' @examples
#' example_celltype_pseudobulks <- generate_synthetic_pseudobulks()
#' conditions <- sapply(strsplit(colnames(example_celltype_pseudobulks), "[_-]"),function(x){x[[2]]})
#' names(conditions) <- colnames(example_celltype_pseudobulks)
#' design <- model.matrix(~conditions)
#' de <- compute_cell_type_specific_DE_parallel(example_celltype_pseudobulks, design, de_method="deseq2")
#' de <- compute_cell_type_specific_DE_parallel(example_celltype_pseudobulks, design, de_method="edger")
#' @export
compute_cell_type_specific_DE_parallel <- function(pseudobulks, design_matrix, fdr=0.05, n.cores=1, de_method=c("deseq2", "edger")) {
	sample <- sapply(strsplit(colnames(pseudobulks), "_"), function(x){x[[1]]})
	cell_type <- sapply(strsplit(colnames(pseudobulks), "_"), function(x){x[[2]]})
	all_outs <- foreach::foreach( type=cell_type) %do% {
		dat <- pseudobulks[,cell_type==type]
		# Create the design matrix for this cell-type
		design <- design_matrix[rownames(design_matrix) %in% colnames(dat),]
		if (Matrix::rankMatrix(design)[1] < ncol(design)) {
			warning(paste("Warning: DE for", type, "could not be computed because the predictors are dependent - design matrix is not full rank."))
		} else {
			if(de_method[1] == "deseq2") {
				list(type=type, de=one_cell_type_DE_edgeR(dat, design, fdr=fdr))
			} else if (de_method[1] == "edger") {
				list(type=type, de=one_cell_type_DE_deseq2(dat, design, fdr=fdr))
			} else {
				print("Error: unrecognized DE method")
			}
		}
	}
	type_names <- sapply(all_outs, function(x){x$type})
	all_outs <- lapply(all_outs, function(x){x$de})
	names(all_outs) <- type_names;
	return(all_outs)
}

#' Cell-type specific DE - one cell-type (edgeR)
#'
#' @description
#' Calculated differential expression for a cell-type given an appropriate design matrix using edgeR
#'
#' @details
#' Called by compute_cell_type_specific_DE to calculate DE of each cell-type. This function filters genes for those with total counts > 5 and detection in at least 3 samples. Then runs DE using edgeR on the provided data & design matrix. Normalization is performed using TMM, and DE is tested using quasi-likelihood (QLF test). Each coefficient is tested and results are combined into a single dataframe. See: edgeR::glmQLFit and glmQLRTest for details.
#' Due to frequent errors in TMM normalization with pseudobulk data, total counts per sample is added to the design matrix automatically.For simplicity, the DE results for total counts are not returned.
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
#' de <- one_cell_type_DE_edgeR(dat, design)
#' @export
one_cell_type_DE_edgeR <- function(dat, design, fdr=0.05) {
	# Filter out lowly expressed genes
	dat <- dat[Matrix::rowSums(dat) > 5 & Matrix::rowSums(dat > 0) >= 3,]
	design <- cbind(design, colSums(dat))
	colnames(design)[ncol(design)] = "totalcount"
	# Perform DE using edgeR
	dge <- edgeR::DGEList(counts=dat)
	dge <- edgeR::calcNormFactors(dge, method="TMM", refColumn=which(Matrix::colSums(dat) == max(Matrix::colSums(dat))))
	dge <- edgeR::estimateDisp(dge, design)
	fit <- edgeR::glmQLFit(dge, design)
	all_contrast_outs <- c()
	for (i in 2:(ncol(design)-1)) {
		contrast_vec <- rep(0, ncol(design))
		contrast_vec[i] <- 1
		de <- edgeR::topTags(edgeR::glmQLFTest(fit, contrast=contrast_vec), n=nrow(dat), p.value=fdr)
		de <- de$table
		de$Coefficient <- colnames(design)[i]
		de$Gene <- rownames(de)
		all_contrast_outs <- rbind(all_contrast_outs, de)
	}
	return(all_contrast_outs)
}

#' Cell-type specific DE - one cell-type (DESeq2)
#'
#' @description
#' Calculated differential expression for a cell-type given an appropriate design matrix using DESeq2
#'
#' @details
#' Called by compute_cell_type_specific_DE to calculate DE of each cell-type. This function filters genes for those with total counts > 5 and detection in at least 3 samples. Then runs DE using deseq2 on the provided data & design matrix. Normalization uses size factors and DE is performed with the default Wald test. see: DESeq2::DESeq for details.
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
#' de <- one_cell_type_DE_deseq2(dat, design)
#' @export
one_cell_type_DE_deseq2 <- function(dat, design, fdr=0.05) {
	# Filter out lowly expressed genes
	dat <- dat[Matrix::rowSums(dat) > 5 & Matrix::rowSums(dat > 0) >= 3,]
	# Perform DE using edgeR
	all_contrast_outs <- c()

	coldata = design;
	dds <- DESeq2::DESeqDataSetFromMatrix(countData = dat,
                              colData = coldata,
                              design=design)
	dds <- DESeq2::DESeq(dds)
	coeffs <- resultsNames(dds) # lists the coefficients
	
	for (i in 2:(length(coeffs))) {
		res <- results(dds, name=coeffs[i])
		# Reformat to look like edgeR results
		de <- res[,c(2,1,4,5,6)]
		colnames(de) <- c("logFC", "baseMean", "stat", "PValue", "FDR")
		de$Coefficient <- colnames(design)[i]
		de$Gene <- rownames(de)
		all_contrast_outs <- rbind(all_contrast_outs, de)
	}
	return(all_contrast_outs)
}
