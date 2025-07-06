#' Plot Enrichment Heatmap
#'
#' @description
#' Summarizes a set of pathway enrichments as a heatmap.
#'
#' @details
#' Selects the top pathways from each set of enrichment results. Then aggregates the significance and enrichment scores for each of the selected pathways across all enrichment results.  
#' @param list_of_rich a list of pathway enrichments for a single DE test obtained from do_ora or do_fgsea. All enrichments should be calculated using either do_ora or do_fgsea.
#' @param pathways a vector of pathway names of pathways to be plot.
#' @param pathways.names a vector of alias names for each pathway to be plot - used to shorted long pathway names when using predetermined pathways, overrides remove.prefix.
#' @param ntop number of pathways to select from each enrichement results.
#' @param colors a vector of colors for the heatmap
#' @param stars a boolean for whether to include stars to indicate significance. 
#' @param stars.col the text colour of the stars
#' @param remove.prefix whether to remove the prefix from each pathway name. 
#' @param prefix.delim the symbol that defines the end of the prefix
#' @param log.scale whether to log (base2) scale the enrichment scores
#' @param bounds boundarys for the colour scheme, if NULL then uses a centered boundaries for the maximum range of values.
#' @param cluster_rows whether to cluster the rows of the heatmap or not
#' @param cluster_cols whether to cluster the columns of the heatmap or not
#' @param plot.result whether to make the plot or just calculate the summarized results.
#' @param both.dir whether to select the ntop pathways both up & down or just the ntop based on pvalue regardless of direction 
#' @param col.anno dataframe of annotations for each comparison. # Not implemented yet
#' @param row.anno dataframe of annotations for each pathway.  # Not implemented yet
#' @return invisibly the data matrix and p value matrix used to make the heatmap (as a list)
#' @examples
#' list_of_rich <- lapply(rpois(5, lambda=10), generate_synthetic_enrichments, ngenes=100)
#' names(list_of_rich) <- paste("celltype", 1:length(list_of_rich), sep="")
#' plot_enrichments_heatmap(list_of_rich)
#' @export

plot_enrichments_heatmap <- function(list_of_rich, pathways=NULL, pathways.names=NULL, ntop=5, colors=grDevices::colorRampPalette(rev(c("red", "orange", "yellow", "black", "navy", "purple", "magenta")))(100), stars=TRUE, stars.col="white", remove.prefix=TRUE, prefix.delim="_", log.scale=FALSE, bounds=NULL, cluster_rows=TRUE, cluster_cols=FALSE, plot.result=TRUE, both.dir=FALSE, col.anno=NULL, row.anno=NULL) {
	# Collect pathways
	all_pathways <- c();
	all_pathways_names <- c();
	if (!is.null(pathways)) {
		all_pathways = pathways
		all_pathway_names = pathways.names
		if (is.null(all_pathway_names)) {
			all_pathway_names = all_pathways
		}
	} else {
		#pathway_origin <- c();
		for (i in names(list_of_rich)) {
			paths <- list_of_rich[[i]]$results
			if (nrow(paths) == 0) {next;}
			paths <- paths[order(paths$FDR),]
			if (both.dir) {
				n_up <- min(nrow(paths[paths[,3] > 0,]), ntop)
				all_pathways <- c(all_pathways, paths[paths[,3] > 0,"pathway"][1:n_up])
				n_dn <- min(nrow(paths[paths[,3] < 0,]), ntop)
				all_pathways <- c(all_pathways, paths[paths[,3] < 0,"pathway"][1:n_dn])

			} else {
				n <- min(nrow(paths), ntop)
				all_pathways <- c(all_pathways, paths$pathway[1:n])
			}
			#pathway_origin <- c(pathway_origin, rep(i, n))
		}
	}
	# remove duplicates
	dups <- duplicated(all_pathways)
	all_pathways <- all_pathways[!dups]
	#pathway_origin <- pathway_origin[!dups]

	# Generate data for heatmap	
	heat_data_all <- generate_heatmap_data(list_of_rich, all_pathways)
	heat_data <- heat_data_all$scores
	heat_pvals <- heat_data_all$pvalues

	# Tidy up pathway names.
	if (plot.result) {
		db=row.anno
		# add annotation for origin of the pathways.
		# anno <- data.frame(origin=pathway_origin)
		# No i don't think this is that useful since there are duplicates that we just remove above...
		if (remove.prefix & is.null(pathways.names)) {
		# Convert removed db tags to a colour label
			pathway_prefixes <- sapply(strsplit(all_pathways, prefix.delim), function(x){x[[1]]})
			all_pathway_names <- sub(paste0("^[^",prefix.delim,"]*",prefix.delim, sep=""), "", all_pathways)
			while(sum(duplicated(all_pathway_names)) > 0) {
				all_pathway_names[duplicated(all_pathway_names)] <- paste0(all_pathway_names[duplicated(all_pathway_names)], "1")
			}
		
			anno_row <- data.frame(pathway_prefixes); 
			rownames(anno_row) <- all_pathway_names;
			if (!is.null(row.anno)) {
				anno_row <- merge(anno_row, row.anno)
			}
			
			rownames(heat_data) <- all_pathway_names
			plot_heatmap(heat_data, heat_pvals, 
					log.scale=log.scale, bounds=bounds, colors=colors, 
					stars=stars, stars.col=stars.col,
					cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
					annotation_row = anno_row, annotation_col = NA)
		} else {
			rownames(heat_data) <- all_pathway_names
			plot_heatmap(heat_data, heat_pvals, 
					log.scale=log.scale, bounds=bounds, colors=colors, 
					stars=stars, stars.col=stars.col,
					cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
					annotation_row = NA, annotation_col = NA)
		}
	}
	invisible(heat_data_all)
}


#' Collect Heatmap Data
#'
#' @description
#' Collects enrichment scores and FDR p-values from a set of enrichments.
#' 
#' @details
#' Combines enrichment scores and FDRs from a list of enrichment output into two matrices which can be plotted as a heatmap across all cell-types.
#' pathways missing from one or another cell-type are assigned scores of 0, and p-values of 1.
#'
#' @param list_of_rich a list of output from do_ora or do_fgsea typically from multiple cell-types
#' @param pathways which pathways to aggregate.
#' @return a list of "scores" and "pvalues" for the provided pathways across all cell-types.
#' @examples
#' list_of_rich <- lapply(rpois(5, lambda=10), generate_synthetic_enrichments, ngenes=100)
#' names(list_of_rich) <- paste("celltype", 1:length(list_of_rich), sep="")
#' heat_data <- generate_heatmap_data(list_of_rich, pathways=c("pathway1", "pathway2", "pathway5", "pathway7"))
#' dim(heat_data$scores)
#' dim(heat_data$pvalues)
#' @export

# list of rich is a list of "out" from our enrichment methods - i.e. contains $results and $contrib
generate_heatmap_data <- function(list_of_rich, pathways) {
	heat_toplot <- matrix(-1, nrow=length(pathways), ncol=length(list_of_rich))
	rownames(heat_toplot) <- pathways
	colnames(heat_toplot) <- names(list_of_rich)
	heat_pvals <- matrix(-1, nrow=length(pathways), ncol=length(list_of_rich))
	rownames(heat_pvals) <- pathways
	colnames(heat_pvals) <- names(list_of_rich)
	for (this_i in names(list_of_rich)) {
		these_rich <- list_of_rich[[this_i]]
		scores <- these_rich$results[match(pathways, these_rich$results$pathway),3]
		pvals <- these_rich$results[match(pathways, these_rich$results$pathway),"FDR"]
		heat_toplot[,this_i] <- scores
		heat_pvals[,this_i] <- pvals
	}
	heat_toplot[is.na(heat_toplot)] <- 0
	heat_pvals[is.na(heat_pvals)] <- 1
	return(list(scores=heat_toplot, pvalues=heat_pvals))
}

#' Plot Heatmap
#'
#' @description
#' Plots a heatmap of a provided enrichment score matrix.
#' 
#' @details
#' log.scale : data is log2-scaled as log2(abs(data)+1)*sign(data) this preserves the +ve/-ve sign which is used for the direction of the expression of the pathway in each cell-type.
#'
#' stars : plots a "*" for enrichments significant at p < 0.05, "**" for p < 0.0005, and "***" for p < 0.000005
#' 
#' @param data a matrix of enrichment scores for various pathways across multiple cell-types
#' @param pvals a matrix of FDR or other adjusted p-values for each enrichment in data
#' @param log.scale whether to log-scale the provided scores (data - See Details).
#' @param bounds the boundaries uses for the color scheme of the enrichment scores. By default will use a centered boundary covering the full range of data.
#' @param colors a vector of colors, the color scale used for enrichment scores.
#' @param stars boolean, whether to plot stars indicating signifcance of the enrichments (see Details).
#' @param stars.col the color for the stars
#' @param cluster_rows whether to cluster & rearrange the rows of the heatmap
#' @param cluster_cols whether to cluster & rearrange the columns of the heatmap
#' @param annotation_row a matrix of annotations for the rows of the heatmaps (see:pheatmap)
#' @param annotation_col a matrix of annotations for the columns of the heatmaps (see:pheatmap)
#' @return Nothing.
#' @examples
#' heat_data <- matrix(rnorm(100), ncol=10)
#' rownames(heat_data) <- paste("pathway", 1:nrow(heat_data), sep="")
#' colnames(heat_data) <- paste("celltype", 1:ncol(heat_data), sep="")
#' pvals=matrix(exp(-abs(rnorm(100,mean=5, sd=10))),ncol=10)
#' plot_heatmap(heat_data, stars=FALSE)
#' plot_heatmap(heat_data, pvals=pvals, stars=TRUE)
#' @export
plot_heatmap <- function(data, pvals=NULL, log.scale=FALSE, bounds=NULL, colors=grDevices::colorRampPalette(rev(c("red", "orange", "yellow", "black", "navy", "purple", "magenta")))(100), stars=TRUE, stars.col="white", cluster_rows=TRUE, cluster_cols=TRUE, annotation_row = NA, annotation_col = NA) {

	if (log.scale) {
		tmp <- sign(data)
		rescaled <- log2(abs(data)+1)*tmp
		data <- rescaled
	}

	if (is.null(bounds)) {
		max_range <- max(abs(data))
		bounds = c(-max_range, max_range)
	}
	star_mat <- matrix("",ncol=ncol(data), nrow=nrow(data))
	if (stars) {
		if(is.null(pvals)) {stop("Error: you must provide an adjusted p-value matrix to plot stars.")}
		if(!identical(dim(data), dim(pvals))) {stop("Error: pvals must be the same dimensions as data.")}
		star_mat[pvals < 0.05] <- "*"
		star_mat[pvals < 0.0005] <- "**"
		star_mat[pvals < 0.000005] <- "***"
	}
	
	pheatmap::pheatmap(data, color=colors, 
		breaks=seq(from=bounds[1], to=bounds[2]+1, length=length(colors)+1),
		display_numbers=star_mat, number_color=stars.col,
		cluster_rows=cluster_rows, cluster_cols=cluster_cols, 
		annotation_row=annotation_row, annotation_col=annotation_col)
}
