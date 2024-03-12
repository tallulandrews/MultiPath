plot_enrichments_heatmap <- function(list_of_rich, ntop=5, colors=colorRampPalette(rev(c("red", "orange", "yellow", "black", "navy", "purple", "magenta")))(100), stars=NULL, stars.col="black", remove.prefix=FALSE, prefix.delim="_", log.scale=FALSE, bounds=NULL, cluster_rows=TRUE, cluster_cols=FALSE, plot.result=TRUE) {
	# Collect pathways
	all_pathways <- c();
	#pathway_origin <- c();
	for (i in names(list_of_rich)) {
		paths <- list_of_rich[[i]]$results
		if (nrow(paths) == 0) {next;}
		paths <- paths[order(paths$FDR),]
		n <- min(nrow(paths), ntop)
		all_pathways <- c(all_pathways, paths$pathway[1:n])
		#pathway_origin <- c(pathway_origin, rep(i, n))
	}
	# remove duplicates
	dups <- duplicated(all_pathways)
	all_pathways <- all_pathways[!dups]
	#pathway_origin <- pathway_origin[!dups]

	# Generate data for heatmap	
	heat_data <- generate_heatmap_data(list_of_rich, all_pathways)

	# Tidy up pathway names.
	if (plot.result) {
		tags=NA
		# add annotation for origin of the pathways.
		# anno <- data.frame(origin=pathway_origin)
		# No i don't think this is that useful since there are duplicates that we just remove above...
		if (remove.prefix) {
		# Convert removed tags to a colour label
			tags <- sapply(strsplit(all_pathways, prefix.delim), function(x){x[[1]]})
			all_pathway_names <- sub(paste0("^[^",prefix.delim,"]*",prefix.delim, sep=""), "", all_pathways)
			anno_row <- data.frame(tags); rownames(anno_row) <- all_pathways;
			rownames(heat_data) <- all_pathway_names
			plot_heatmap(heat_data, log.scale=log.scale, bounds=bounds, colors=colors, 
					cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
					annotation_row = anno_row, annotation_col = NA)
		} else {
			all_pathway_names <- all_pathways
			rownames(heat_data) <- all_pathway_names
			plot_heatmap(heat_data, log.scale=log.scale, bounds=bounds, colors=colors, 
					cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
					annotation_row = NA, annotation_col = NA)
		}
	}
	return(heat_data)
}

# list of rich is a list of "out" from our enrichment methods - i.e. contains $results and $contrib
generate_heatmap_data <- function(list_of_rich, pathways) {
	heat_toplot <- matrix(-1, nrow=length(pathways), ncol=length(list_of_rich))
	rownames(heat_toplot) <- pathways
	colnames(heat_toplot) <- names(list_of_rich)
	for (this_i in names(list_of_rich)) {
		these_rich <- list_of_rich[[this_i]]
		scores <- these_rich$results[match(pathways, these_rich$results$pathway),3]
		heat_toplot[,this_i] <- scores
	}
	heat_toplot[is.na(heat_toplot)] <- 0
	return(heat_toplot)
}

plot_heatmap <- function(data, log.scale=FALSE, bounds=NULL, colors=colorRampPalette(rev(c("red", "orange", "yellow", "black", "navy", "purple", "magenta")))(100), cluster_rows=TRUE, cluster_cols=TRUE, annotation_row = NA, annotation_col = NA) {

	if (log.scale) {
		tmp <- sign(data)
		rescaled <- log2(abs(data)+1)*tmp
		data <- rescaled
	}

	if (is.null(bounds)) {
		max_range <- max(abs(data))
		bounds = c(-max_range, max_range)
	}

	pheatmap::pheatmap(data, color=colors, 
		breaks=seq(from=bounds[1], to=bounds[2], length=length(colors)+1),
		cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
		annotation_row=annotation_row, annotation_col=annotation_col)
}
