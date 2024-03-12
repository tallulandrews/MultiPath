# I have only tested this after running the following:
# library(msigdb)
# library(ExperimentHub)
# library(GSEABase)

get_pathways <- function(species=c("mouse", "human"), include.celltype=FALSE) {
	if (species == "mouse") {
		out <- msigdb::getMsigdb(org="mm", id="SYM")
		if (include.celltype) {
			genesets <- subsetCollection(out, c('c2','c5','h', 'c8'), c('GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
		} else {
			genesets <- subsetCollection(out, c('c2','c5','h'), c('GO:BP', 'CP:REACTOME', 'CP:BIOCARTA'))
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
	if (exists("genesets")) {
		return(genesets)
	} else {
		print(paste("Did not recognize:", species, "no genesets found."))
	}
}

Convert_GSEAObj_to_List <- function(gsea_obj) {
	out <- NULL
	if (class(gsea_obj)[1] =="GeneSet") {
		out <- list(gsea_obj@geneIds)
		names(out) <- gsea_obj@setName
		return(out)
	}
	if(class(gsea_obj)[1] == "GeneSetCollection") {
		tmp <- lapply(gsea_obj, Convert_GSEAObj_to_List)
		return(sapply(tmp, function(x){x}))
	}
}
