# I have only tested this after running the following:
# library(msigdb)
# library(ExperimentHub)
# library(GSEABase)

#' Get Pathways from MSigDb
#'
#' @description
#' Downloads common pathways from msigdb for either human or mouse. 
#' 
#' @details
#' Downloads MSigDB Hallmark pathways, GO biological processes, Reactome and Biocarta pathways from `msigdb` for either human or mouse. For human it also adds the KEGG pathways. Optionally also downloads cell-type marker genesets.
#' @param species whether to download human or mouse pathways
#' @param include.celltype whether to include cell-type marker genesets.
#' @return GeneSetCollection object of pathways.
#' @examples
#' genesets <- get_pathways("mouse")
#' genesets <- get_pathways("human")
#' @export
get_pathways <- function(species=c("mouse", "human"), include.celltype=FALSE) {
	if (species == "mouse") {
		out <- msigdb::getMsigdb(org="mm", id="SYM")
		if (include.celltype) {
			genesets <- subsetCollection(out, c('h', 'c8'), c('GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
		} else {
			genesets <- subsetCollection(out, c('h'), c('GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
		}
	}
	if (species == "human") {
		out <- msigdb::getMsigdb(org="hs", id="SYM")
		out <- appendKEGG(out)
		# collections: c8 = cell type signatures
		#		h = hallmark gene sets
		#		c5 = GO
		#		c2 = curated pathways
		if (include.celltype) {
			genesets <- subsetCollection(out, collection=c('h','c8'), subcollection=c('CP:KEGG', 'GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
		} else {
			genesets <- subsetCollection(out, collection=c('h'), subcollection=c('CP:KEGG', 'GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
		}
	}
	if (species == "test") { # for test-case for this package
		out <- msigdb::getMsigdb(org="hs", id="SYM")
		genesets <- subsetCollection(out, collection=c('h'))
	}
	if (exists("genesets")) {
		return(genesets)
	} else {
		print(paste("Did not recognize:", species, "no genesets found."))
	}
}

#' Convert GSEAObj
#'
#' @description
#' Converts the GSEA object obtained from msigdb or other Bioconductor pathway database to a list suitable for use in `do_ora`.
#' 
#' @param gsea_obj a GeneSet or GeneSetCollection object, such as those obtained from the msigdb database.
#' @return a list of gene sets.
#' @examples
#' genesets <- get_pathways("test")
#' geneset_list <- convert_GSEAObj_to_list(genesets)
#' length(geneset_list)
#' @export
convert_GSEAObj_to_list <- function(gsea_obj) {
	out <- NULL
	if (class(gsea_obj)[1] =="GeneSet") {
		out <- list(gsea_obj@geneIds)
		names(out) <- gsea_obj@setName
		return(out)
	}
	if(class(gsea_obj)[1] == "GeneSetCollection") {
		tmp <- lapply(gsea_obj, convert_GSEAObj_to_list)
		return(sapply(tmp, function(x){x}))
	}
}
