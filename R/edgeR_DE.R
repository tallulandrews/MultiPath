compute_cell_type_specific_DE <- function(pseudobulks, design_matrix, fdr=0.05) {
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

compute_cell_type_specific_DE_parallel <- function(pseudobulks, design_matrix, fdr=0.05) {
	#### THIS IS NOT SET UP YET ####
	sample <- sapply(strsplit(colnames(pseudobulk_expr_mat_all), "_"), function(x){x[[1]]})
	cell_type <- sapply(strsplit(colnames(pseudobulk_expr_mat_all), "_"), function(x){x[[2]]})
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
