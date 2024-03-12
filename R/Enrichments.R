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


# Optimized ORA
# adpated from: https://jokergoo.github.io/2023/04/05/speed-up-over-representation-enrichment-analysis/
# originally created by: Zuguang Gu
	
library(Rcpp)
sourceCpp(code = '
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <unordered_set>
using namespace Rcpp;

// [[Rcpp::export]]
List intersectToList(List lt, StringVector x) {

    int n = lt.size();
    List out(n);

    std::unordered_set<String> seen;
    seen.insert(x.begin(), x.end());

    for(int i = 0; i < n; i++) {
      
        StringVector v = as<StringVector>(lt[i]);
        LogicalVector l(v.size());

        std::unordered_set<String> seen2;

        for(int j = 0; j < v.size(); j ++) {
            l[j] = seen.find(v[j]) != seen.end() && seen2.insert(v[j]).second;
        }

        out[i] = v[l];
    }

    return out;
}
')

do_ora <- function(sig_genes, pathways, background, fdr=0.05, min.term.size=15, max.term.size=1000, seed=2910){
	set.seed(seed)
	# Filter Pathways
	path_size <- sapply(pathways, length)
	pathways <- pathways[path_size < max.term.size & path_size > min.term.size]

	# Run hypergeometric test
	sig_genes <- intersect(sig_genes, background) # ensure only considering genes that exist in the background
	path_names <- names(pathways)
	sig_genes <- intersect(sig_genes, background)
	pathways <- intersectToList(pathways, background) 
	names(pathways) <- path_names

	n_background <- length(background)
	n_genes <- length(sig_genes)
    
	x <- sapply(intersectToList(pathways, sig_genes), length)  # this line has been improved
	m <- sapply(pathways, length)
	n <- n_background - m
	k <- n_genes
    
	p <- phyper(x - 1, m, n, k, lower.tail = FALSE)
	names(p) <- path_names
	fdr_res <- p.adjust(p, method="fdr")

	# Generate nice results
	keep <- fdr_res < fdr
	if (sum(keep) == 0) {warning("Warning:No significant enrichments"); return();}
	pathways <- pathways[keep]
	contrib_genes <- intersectToList(pathways, sig_genes)
	names(contrib_genes) <- names(pathways)
	score <- log2((x[keep]/k)/(m[keep]/n_background))
	res = data.frame(pathway=names(pathways), intersection=x[keep], log2fe=score, FDR=fdr_res[keep]) 
	return(list(results=res, contrib=contrib_genes))
}

# ================== Functions for condensing synonymous pathways ================ #
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


cluster_overlaps <- function(overlaps, equivalent=0.5, plot.result=TRUE) {
	groups <- cutree(hclust(as.dist(1-overlaps), method="single"), h=1-equivalent)
	if (plot.result) {
		unique_paths = names(table(groups))[table(groups) == 1]
		anno <- data.frame(group=groups)
		anno[anno[,1] %in% unique_paths,1] = "unique"
		pheatmap::pheatmap(overlaps, annotation_col=anno, annotation_legend=FALSE, main="Gene Overlaps")	
	}
	return(groups)
}

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

# select one representative for each group of overlapping terms
condense_terms <- function(out, equivalent=0.5, verbose=FALSE, prioritize.signaling=TRUE) {
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

trim_pathway_names <- function(out, prefix_length=1, nwords=4, split_char="_") {
	curr_names <- out$results$pathway
	split <- strsplit(curr_names, split_char)
	out$results$trimmed_name <- sapply(split, function(x){paste(x[(1+prefix_length):min(length(x),nwords+prefix_length)], collapse="_")})
	return(out)
}

