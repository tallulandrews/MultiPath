#' Teat Pathways
#'
#' @description MSigDb Hallmark pathways
#'
#' @format a list of 50 pathways:
#' \describe{
#'    \item{item name}{name of the pathway}
#'    \item{value}{vector of human gene IDs}
#'  }
#'
#' @source {msigdb} Bioconductor database package containing MSigDB pathways.
"test_paths"


#' GeneSet Enrichment Analysis
#'
#' @description
#' Performs gene set enrichment analysis using fgsea
#' 
#' @details
#' Uses Rcpp backend to optimize overrepresentation analysis (ORA). 
#' A given set of query genes is tested for enrichment in genes from provided named pathways or genesets, compared to a specified set of background genes. 
#' Typically a query gene set will be a group of up or down regulated genes obtained from differential expression testing. The background would be all the genes in the count matrix used in the differential expression test. 
#' ORA is performed with a hypergeometric test against a specified background and applies a FDR multiple testing correction, after filtering provided pathways for terms with a reasonable number of genes attributed to it.
#' @param scored_genes a vector of scores named with the gene name.
#' @param pathways a list of gene sets or pathways to test for enrichment. The list must be named by the pathway name and contain vectors of gene names for the genes in that gene set.
#' @param fdr the FDR threshold used to filter output.
#' @param min.term.size pathways must have at least this number of genes to be tested.
#' @param max.term.size pathways must have fewer than this number of genes to be tested.
#' @param seed random seed that will be set at the start of this function, this ensure complete reproducibility.
#' @param nperm the number of permutations used to estimate p-values, it is recommended to use at least 30,000 to ensure reproducible results.
#' @return A list containing two items:
#' results = a dataframe containing the pathway name, number of query genes in the pathway (intersection), log2foldenrichment, FDR
#' contrib = a list of genes found in the intersection between this pathway and the query gene set. 
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths
#' l2fc <- rnorm(100, sd=2)
#' names(l2fc) <- sample(unlist(paths), size=100)
#' rich <- do_gsea(l2fc, paths, nperm=50) # nperm should be a minimum of 50000 for reproducible results
#' @export
# note: nperm is set to 100000 because of random instability seen at nperm = 10000
do_gsea <- function(scored_genes, pathways, fdr=0.05, min.term.size=15, max.term.size=1000, seed=2910, nperm=100000){
	set.seed(seed)
	scored_genes[scored_genes > 0 & !is.finite(scored_genes)] <- max(scored_genes[is.finite(scored_genes)])+1
	scored_genes[scored_genes < 0 & !is.finite(scored_genes)] <- min(scored_genes[is.finite(scored_genes)])-1
	rich <- fgsea::fgseaMultilevel(pathways, scored_genes, minSize=min.term.size, maxSize=max.term.size, nPermSimple=nperm)
	keep <- !is.na(rich$pval) & rich$padj < fdr
	if (sum(keep) == 0) {warning("Warning:No significant enrichments"); return();}
	
	#contrib_genes <- as.list(rich[,8][[1]]) # This doesn't work when in a package????
	contrib_genes <- lapply(1:nrow(rich), function(x){unlist(rich[x,8])})
	names(contrib_genes) <- rich$pathway
	n_contrib <- sapply(contrib_genes, length)
	
	#print(length(rich$pathway))
	#print(length(n_contrib))
	#print(length(rich$NES))
	#print(length(rich$padj))
	#return(rich)
	res = data.frame(pathway=rich$pathway, intersection=n_contrib, NES=rich$NES, FDR=rich$padj)
	return(list(results=res[keep,], contrib=contrib_genes[keep]))
}

#' Optimized Overrepresentation Analysis
#'
#' @description
#' Performs overrepresentation analysis using a hypergeometic test with an FDR correction.
#' 
#' @details
#' Uses Rcpp backend to optimize overrepresentation analysis (ORA). 
#' A given set of query genes is tested for enrichment in genes from provided named pathways or genesets, compared to a specified set of background genes. 
#' Typically a query gene set will be a group of up or down regulated genes obtained from differential expression testing. The background would be all the genes in the count matrix used in the differential expression test. 
#' ORA is performed with a hypergeometric test against a specified background and applies a FDR multiple testing correction, after filtering provided pathways for terms with a reasonable number of genes attributed to it.
#'
#' adpated from: https://jokergoo.github.io/2023/04/05/speed-up-over-representation-enrichment-analysis/
#' by: Zuguang Gu
#' @param sig_genes a vector of gene names that is the query for the ORA
#' @param pathways a list of gene sets or pathways to test for enrichment. The list must be named by the pathway name and contain vectors of gene names for the genes in that gene set.
#' @param background a vector of gene games to be used as the background.
#' @param fdr the FDR threshold used to filter output.
#' @param min.term.size pathways must have at least this number of genes in the background list of genes to be tested.
#' @param max.term.size pathways must have fewer than this number of genes in the background list of genes to be tested.
#' @param include.underrepresented whether to include pathways underrepresented among the provided genes
#' @return A list containing two items:
#' results = a dataframe containing the pathway name, number of query genes in the pathway (intersection), log2foldenrichment, FDR
#' contrib = a list of genes found in the intersection between this pathway and the query gene set. 
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' dim(rich$results)[1] == length(rich$contrib) # TRUE
#' @export
do_ora <- function(sig_genes, pathways, background, fdr=0.05, min.term.size=10, max.term.size=1000, include.underrepresented=FALSE){
	# Filter Pathways
	tmp <- names(pathways)
	pathways <- intersectToList(pathways, background) 
	names(pathways) <- tmp
	path_size <- sapply(pathways, length)
	if (max(path_size) == 0) {stop("Error:No overlap between pathways and background. Please check your genes are stored as appropriate gene symbols.")}
	pathways <- pathways[path_size < max.term.size & path_size > min.term.size]

	# Run hypergeometric test
	sig_genes <- intersect(sig_genes, background) # ensure only considering genes that exist in the background
	path_names <- names(pathways)
	names(pathways) <- path_names

	n_background <- length(background)
	n_genes <- length(sig_genes)
    
	x <- sapply(intersectToList(pathways, sig_genes), length)  # this line has been improved
	m <- sapply(pathways, length)
	n <- n_background - m
	k <- n_genes
    
	p <- stats::phyper(x - 1, m, n, k, lower.tail = FALSE)
	names(p) <- path_names
	fdr_res <- stats::p.adjust(p, method="fdr")

	# Generate nice results
	if (include.underrepresented) {
		keep <- fdr_res < fdr
	} else {
		keep <- fdr_res < fdr & ( (x/k) > (m/n_background) )
	}
	if (sum(keep) == 0) {warning("Warning:No significant enrichments"); return();}
	pathways <- pathways[keep]
	contrib_genes <- intersectToList(pathways, sig_genes)
	names(contrib_genes) <- names(pathways)
	x[x==0] <- 1/n_background
	score <- log2(((x[keep])/k)/(m[keep]/n_background)) # Add the 0.1 to avoid -Inf
	if (n_genes == 0) {
		# Catch cases where no genes provided
		score <- rep(0, length(score))
		fdr_res <- rep(1, length(fdr_res))
		warning("Warning: no genes in query set")
	}
	res = data.frame(pathway=names(pathways), intersection=x[keep], log2fe=score, FDR=fdr_res[keep]) 
	return(list(results=res, contrib=contrib_genes))
}

# ================== Functions for condensing synonymous pathways ================ #

#' Pathway Overlaps
#'
#' @description
#' Calculates the degree of overlap between genes contributing to each pathway enrichment.
#' 
#' @details
#' Overlaps are calculated as the length of the intersection divided by the total genes in the larger pathway.
#' "Max" = number of genes in the largest pathway. <- default
#' "Min" = number of genes in the smallest pathway.
#' "Union" = number of genes in the union of the two pathways.
#' @param out a set of pathway enrichments for a single DE test obtained from do_ora or do_fgsea
#' @param denom one of "max", "min" or "avg" determining how the denominator is calculated (See: details)
#' @return A matrix of overlap scores for all pairs of pathways.
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths;print("1")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)));
#' overlaps <- get_overlaps(rich)
#' @export
get_overlaps <- function(out, denom=c("max", "min", "union")) {
	total_genes <- sapply(out$contrib, length)
#	olap <- lapply(out$contrib, function(x) {sapply(intersectToList(out$contrib, x), length)})
#	olap <- matrix(unlist(olap), ncol=length(out$contrib))
#	rownames(olap) <- colnames(olap) <- names(out$contrib)

	# Rcpp-ized
	olap <- overlapsListvsList(out$contrib, out$contrib)
	rownames(olap) <- colnames(olap) <- names(out$contrib)

	# get denominator - size of the largest pathway
	tab <- cbind(rep(total_genes, length(total_genes)), rep(total_genes, each=length(total_genes)))
	if (denom[1] == "max") {
		denom <- matrix(apply(tab, 1, max), ncol=length(total_genes))
	} else if (denom[1] == "min") {
		denom <- matrix(apply(tab, 1, min), ncol=length(total_genes))
	} else if (denom[1] == "union") {
		denom <- matrix(apply(tab, 1, max), ncol=length(total_genes)) + 
			 matrix(apply(tab, 1, min), ncol=length(total_genes)) -
			 olap
	}

	return(olap/denom)
}


#' Cluster Pathways
#'
#' @description
#' Clusters pathways based on overlapping contributing genes.
#' 
#' @details
#' Uses single-linkage hierarchical clustering(`hclust`) to group pathway terms based on overlapping contributing genes.
#' @param overlaps a matrix of pair-wise overlaps between contributing genes to each pathway, this is obtained from 'get_overlaps'.  
#' @param equivalent a scalar value for the equivalence threshold.
#' @param plot.result a boolean determining whether to plot a heatmaps of the overlaps showing he pathway groups.
#' @return A vector of group IDs for each pathway.
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths;print("2")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)));
#' overlaps <- get_overlaps(rich);
#' pathway_groups <- cluster_overlaps(overlaps)
#' @export
cluster_overlaps <- function(overlaps, equivalent=0.5, plot.result=TRUE) {
	groups <- stats::cutree(stats::hclust(stats::as.dist(1-overlaps), method="single"), h=1-equivalent)
	names(groups) <- colnames(overlaps)
	if (plot.result) {
		unique_paths = names(table(groups))[table(groups) == 1]
		anno <- data.frame(group=groups)
		anno[anno[,1] %in% unique_paths,1] = "unique"
		pheatmap::pheatmap(overlaps, annotation_col=anno, annotation_legend=FALSE, main="Gene Overlaps")	
	}
	return(groups)
}

#' Select Term
#'
#' @description
#' Selects one pathway to represent a group of pathways
#' 
#' @details
#' Given a set of synonymous pathways, this function selects one pathway to represent the group. Each pathway is scored based on the length of the pathway name, the p-value of the pathway enrichment, the size of the intersection between the pathway and the DE gens, and whether or not the pathway is a signalling pathway. Each value is ranked across pathways using the "min" method to resolve ties.
#' These ranks are summed and the best scoring pathway is selected.
#' @param out a set of pathway enrichments for a single DE test obtained from do_ora or do_fgsea
#' @param terms a vector of names of pathways that are synonymous
#' @param path_scores a dataframe containing the columns: "pathway", "Score", "FDR". This is used primarily used as part of condensed_terms_multi, to prioritize pathways first by the Score, then in the case of tied Scores by FDR.  
#' @param verbose boolean, whether to print data table used as the basis for selecting the term.
#' @return The name of one pathway to represent the group.
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths;print("3")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)));
#' overlaps <- get_overlaps(rich);
#' pathway_groups <- cluster_overlaps(overlaps)
#' chosen1 <- select_term(rich, rich$results[pathway_groups == "1","pathway"])
#' chosen2 <- select_term(rich, rich$results[,"pathway"])
#' @export
select_term <- function(out, terms, verbose=FALSE, path_scores=NULL) {
	#term_length <- sapply(strsplit(terms, "[ _]"), length)
	#intersection_size <- out$results$intersection[match(terms, out$results$pathway)]
	pvalue <- out$results$FDR[match(terms, out$results$pathway)]
	#is_signalling <- grepl("signaling", terms, ignore.case=TRUE)+0

	# Key words: Signaling
	# intersection size : larget == better
	# length of term name : smaller == better
	# pvalue : smaller == better
	# consensus <- rank(term_length, ties.method="min")+rank(1/intersection_size, ties.method="min")+rank(pvalue, ties.method="min")+rank(1-is_signalling, ties.method="min")
	# chosen <- consensus == min(consensus)
	chosen <- pvalue == min(pvalue)
	if (!is.null(path_scores)) {
		if (! sum(colnames(path_scores) %in% c("pathway", "Score", "FDR"))==3) {
			stop("Error: path_scores does not have the required columns.") 
		}
		path_scores2 <- path_scores[match(terms, path_scores[,"pathway"]),]
		if (!identical(path_scores2[,"pathway"],terms)) {
			stop('Error: The "pathway" column of the path_scores table does not match the provided pathway enrichments.')
		}
		# First choose by maximum Score
		max_score <- max(path_scores2[,"Score"])
		top_candidates <- path_scores2[,"Score"] == max_score
		if (sum(top_candidates) > 1) {
			# Then choose by lowest avg FDR
			min_FDR <- min(path_scores2[top_candidates,"FDR"])
			chosen <- top_candidates & path_scores2[,"FDR"] <= min_FDR
		} else {
			chosen <- top_candidates
		}
	}
	if (verbose){
		gene_frequencies <- table(unlist(out$contrib[names(out$contrib) %in% terms]))
		#print(data.frame(terms, term_length, intersection_size, pvalue, is_signalling, consensus, chosen))
		print(data.frame(terms, pvalue, chosen))
		print(gene_frequencies)
	}
	if (sum(chosen) == 0) {
		stop("Error: Failed to chose a representative term")
	}
	return(terms[chosen][1]) # ensure only one term is returned for each group.
}

#' Condense Terms
#'
#' @description
#' Removes synonymous pathways from a set of output pathway enrichments
#' 
#' @details 
#' First calculates overlaps using 'get_overlaps' then identifies groups of synonymous pathways using 'cluster_overlaps' then selects the best representative term for each group using 'select_term'.
#'
#' @param out a set of pathway enrichments for a single DE test obtained from do_ora or do_fgsea
#' @param equivalent a scalar value for the equivalence threshold.
#' @param verbose boolean, whether to print data table used as the basis for selecting representative terms.
#' @param path_scores a dataframe containing the columns: "pathway", "Score", "FDR". This is used primarily used as part of condensed_terms_multi, to prioritize pathways first by the Score, then in the case of tied Scores by FDR.  
#' @return The same structured data like the output from do_ora or do_fgse but with synonmyous termscondensed into a single pathway
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths;print("4")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' rich <- condense_terms(rich)
#' @export

# select one representative for each group of overlapping terms
condense_terms <- function(out, equivalent=0.5, verbose=FALSE, path_scores=NULL) {
	if (is.null(out)) {
		warning("Warning: No pathways provided to condense_terms.")
		return(out);	
	}
	if (is.null(dim(out$result))) {
		return(out)
	}
	if(nrow(out$result) == 1) {
		return(out)
	}
	overlaps <- get_overlaps(out)
	groups <- cluster_overlaps(overlaps, equivalent=equivalent, plot.result=verbose)
	non_unique <- names(table(groups))[table(groups)>1]
	keep <- names(groups)[!groups %in% non_unique]
	for (g in non_unique) {
		chosen_term <- select_term(out, names(groups)[groups==g], verbose=verbose, path_scores=path_scores)
		keep <- c(keep, chosen_term)
	}
	new_out <- list(results=out$results[out$results$pathway %in% keep,], contrib=out$contrib[names(out$contrib) %in% keep])
	return(new_out)
}

#' Condense Terms Multi
#'
#' @description
#' Removes synonymous pathways from multiple sets of output pathway enrichments
#' 
#' @details 
#' Combines multiple sets of pathways enrichments for different cell-types / conditions 
#' that have not yet been condensed. Then finds overlaps between all pathway enrichments. 
#' Finally identifies a since consensus pathway to represent each cluster of synonymous terms.
#' 
#' Consensus pathway is defined as: ....
#' See: \code{cluster_overlaps} for details of pathway clustering.
#' See: \code{get_oberlaps} for deatils of overlap calculation.
#'
#' @param out_list a list of pathway enrichments for multiple DE tests obtained from do_ora or do_fgsea
#' @param equivalent a scalar value for the equivalence threshold. (see: \code{get_overlaps})
#' @param verbose boolean, whether to print data table used as the basis for selecting representative terms.
#' @return The same structured list of pathway enrichments but with redundant terms removed.
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths;print("5")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)));
#' rich <- condense_terms(rich)
#' @export

# select one representative for each group of overlapping terms
# This is super memory inefficient
# Use Rowan Canario's New Approach Instead
condense_terms_multi <- function(out_list, equivalent=0.5, verbose=FALSE) {

	# Calculate frequency and average FDR for each pathway across all tests
	all_pathways <- lapply(out_list, function(x) x[["results"]][,"pathway"])
	pathway_freq <- table(unlist(all_pathways))
	all_FDRs <- lapply(out_list, function(x) x[["results"]][,"FDR"])
	avg_FDRs <- stats::aggregate(unlist(all_FDRs), by=list(unlist(all_pathways)), mean)
	if (!identical(avg_FDRs[,1], names(pathway_freq))) {
		stop("Error: FDRs do not match Freqs in condense_terms_multi")
	}
	pathway_dat <- data.frame(pathway=names(pathway_freq), Score=as.numeric(unlist(pathway_freq)), FDR=avg_FDRs[,2])

	for (i in 1:length(out_list)) {
		out_list[[i]] <- condense_terms(out_list[[i]], path_scores=pathway_dat)
	}
	return(out_list)
}



obsolete_condense_terms_multi <- function(out, equivalent=0.5, verbose=FALSE) {
	separated_symbol = ";="
	if (is.null(out)) {
		warning("Warning: No pathways provided to condense_terms_multi.")
		return(out);	
	}
	if(length(out) == 1) {
		return(condense_terms(out, equivalent=equivalent, verbose=verbose))
	}
	# Remove any thing in out that has no enrichments
	check <- sapply(out,function(x){dim(x$results)})
	out <- out[!is.null(check)] # NOTE: This also removes those with only 1 enrichment

	# merge all the enrichments into a single thing
	merged_enrichments <- list(results=c(), contrib=list())
	for (i in 1:length(out)) {
		this_name <- names(out)[i]
		out[[i]]$results$condition <- this_name
		merged_enrichments$results <- rbind(merged_enrichments$results, out[[i]]$results)
		names(out[[i]]$contrib) <- paste(this_name, names(out[[i]]$contrib), sep=separated_symbol)
		merged_enrichments$contrib <- append( merged_enrichments$contrib, out[[i]]$contrib)
	}
	# Identify equivalent terms based on clustering of overlaps between them
	overlaps <- get_overlaps(merged_enrichments)
	groups <- cluster_overlaps(overlaps, equivalent=equivalent, plot.result=verbose)

	# Identify unique terms, these we keep as is.
	non_unique <- names(table(groups))[table(groups)>1]
	keep <- names(groups)[!groups %in% non_unique]
	tmp <- strsplit(keep, separated_symbol)
	keep <- sapply(tmp, function(x){x[2]})	
	# For each non-unique term, we need to identify a consensus term
	# Criteria: 
	# - Should be found in all of the conditions that the group of pathways is found in
	# - Should be the most significant on average across all the conditions.
	for (g in non_unique) {
		synonymous_terms <- names(groups)[groups==g]
		stuff <- strsplit(synonymous_terms, separated_symbol) 
		terms <- sapply(stuff, function(x){x[2]})
		condition <- sapply(stuff, function(x){x[1]})
		# Number of conditions each term is found in
		n_condition <- stats::aggregate(condition, by=list(terms), function(x){length(unique(x))})
		potential <- n_condition[n_condition[,2]==max(n_condition[,2]),1]
		
                # average p-values
		stats <- merged_enrichments$results[merged_enrichments$results$pathway %in% potential,]
		overall_score <- stats::aggregate(stats$FDR, by=list(stats$pathway), mean)
		chosen_term <- overall_score[overall_score$x == min(overall_score$x),1]
		keep <- c(keep, chosen_term)
	}
	new_out <- list()
	for (i in 1:length(out)) {
		this_name <- names(out)[i]
		new_out[[this_name]] <- list(results=out[[i]]$results[out[[i]]$results$pathway %in% keep,], 
					     contrib=out[[i]]$contrib[names(out[[i]]$contrib) %in% keep])
		res <- out[[i]]$results
		contrib <- out[[i]]$contrib
		names(contrib) <- sapply(strsplit(names(contrib), separated_symbol), function(x){x[[2]]})

		new_out[[this_name]] <- list(results=res[res$pathway %in% keep,], 
					     contrib=contrib[names(contrib) %in% keep])
		if (nrow(new_out[[this_name]]$results) != length(new_out[[this_name]]$contrib)) {warning("What? Something's wrong.")}
	}
	return(new_out)
}


#' Trim Pathway Names
#'
#' @description
#' Trims the length of pathway names.
#' 
#' @details 
#' Removes a prefix and truncates the pathway names to a specified maximum number of words. 
#' This improves plot appearance when very long pathways are selected.
#'
#' @param out a set of pathway enrichments for a single DE test obtained from do_ora or do_fgsea
#' @param prefix_length the number of 'words' that make up the prefix you want to remove.
#' @param nwords the maximum number of words in the final pathway name.
#' @param split_char the regular expression to use to split the pathway name into individual words
#' @return pathway enrichments provided in the "out" argument with an additional column in the results that includes the trimmed pathway name.
#' @examples
#' #paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' paths <- test_paths; print("6")
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' rich <- condense_terms(rich)
#' rich <- trim_pathway_names(rich)
#' @export
trim_pathway_names <- function(out, prefix_length=1, nwords=4, split_char="_") {
	curr_names <- out$results$pathway
	split <- strsplit(curr_names, split_char)
	out$results$trimmed_name <- sapply(split, function(x){paste(x[(1+prefix_length):min(length(x),nwords+prefix_length)], collapse="_")})
	return(out)
}

