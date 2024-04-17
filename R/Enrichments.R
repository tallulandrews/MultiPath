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
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
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
	
	contrib_genes <- as.list(rich[,8][[1]])
	n_contrib <- sapply(contrib_genes, length)
	
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
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' dim(rich$results)[1] == length(rich$contrib) # TRUE
#' @export
do_ora <- function(sig_genes, pathways, background, fdr=0.05, min.term.size=10, max.term.size=1000, include.underrepresented=FALSE){
	# Filter Pathways
	tmp <- names(pathways)
	pathways <- intersectToList(pathways, background) 
	names(pathways) <- tmp
	path_size <- sapply(pathways, length)
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
	score <- log2(((x[keep]+0.1)/k)/(m[keep]/n_background)) # Add the 0.1 to avoid -Inf
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
#' @param out a set of pathway enrichments for a single DE test obtained from do_ora or do_fgsea
#' @return A matrix of overlap scores for all pairs of pathways.
#' @examples
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' overlaps <- get_overlaps(rich)
#' @export
get_overlaps <- function(out) {
	total_genes <- sapply(out$contrib, length)
	olap <- lapply(out$contrib, function(x) {sapply(intersectToList(out$contrib, x), length)})
	olap <- matrix(unlist(olap), ncol=length(out$contrib))
	rownames(olap) <- colnames(olap) <- names(out$contrib)

	# get denominator - size of the largest pathway
	tab <- cbind(rep(total_genes, length(total_genes)), rep(total_genes, each=length(total_genes)))
	denom <- matrix(apply(tab, 1, max), ncol=length(total_genes))

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
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' overlaps <- get_overlaps(rich)
#' pathway_groups <- cluster_overlaps(overlaps)
#' @export
cluster_overlaps <- function(overlaps, equivalent=0.5, plot.result=TRUE) {
	groups <- stats::cutree(stats::hclust(stats::as.dist(1-overlaps), method="single"), h=1-equivalent)
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
#' @param verbose boolean, whether to print data table used as the basis for selecting the term.
#' @param prioritize.signaling a boolean, indicating whether pathways recieve a bonus for being signalling pathways
#' @return The name of one pathway to represent the group.
#' @examples
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' overlaps <- get_overlaps(rich)
#' pathway_groups <- cluster_overlaps(overlaps)
#' chosen1 <- select_term(rich, rich$pathway[pathway_groups == "1"])
#' @export
select_term <- function(out, terms, verbose=FALSE, prioritize.signaling=TRUE) {
	term_length <- sapply(strsplit(terms, "[ _]"), length)
	intersection_size <- out$results$intersection[match(terms, out$results$pathway)]
	pvalue <- out$results$FDR[match(terms, out$results$pathway)]
	is_signalling <- grepl("signaling", terms, ignore.case=TRUE)+0

	# Key words: Signaling
	# intersection size : larget == better
	# length of term name : smaller == better
	# pvalue : smaller == better
	consensus <- rank(term_length, ties.method="min")+rank(1/intersection_size, ties.method="min")+rank(pvalue, ties.method="min")+rank(1-is_signalling, ties.method="min")
	chosen <- consensus == min(consensus)
	if (verbose){
		gene_frequencies <- table(unlist(out$contrib[names(out$contrib) %in% terms]))
		print(data.frame(terms, term_length, intersection_size, pvalue, is_signalling, consensus, chosen))
		print(gene_frequencies)
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
#' @param prioritize.signaling a boolean, indicating whether pathways recieve a bonus for being signalling pathways
#' @return The same structured data like the output from do_ora or do_fgse but with synonmyous termscondensed into a single pathway
#' @examples
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
#' rich <- do_ora(paths[[1]], paths, background=unique(unlist(paths)))
#' rich <- condense_terms(rich)
#' @export

# select one representative for each group of overlapping terms
condense_terms <- function(out, equivalent=0.5, verbose=FALSE, prioritize.signaling=TRUE) {
	if (is.null(out)) {
		warning("Warning: No pathways provided to condense_terms.")
		return(out);	
	}
	if(nrow(out$result) == 1) {
		return(out)
	}
	overlaps <- get_overlaps(out)
	groups <- cluster_overlaps(overlaps, equivalent=equivalent, plot.result=verbose)
	non_unique <- names(table(groups))[table(groups)>1]
	keep <- names(groups)[!groups %in% non_unique]
	for (g in non_unique) {
		chosen_term <- select_term(out, names(groups)[groups==g], verbose=verbose, prioritize.signaling=prioritize.signaling)
		keep <- c(keep, chosen_term)
	}
	new_out <- list(results=out$results[out$results$pathway %in% keep,], contrib=out$contrib[names(out$contrib) %in% keep])
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
#' paths <- convert_GSEAObj_to_list(get_pathways("test"))
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

