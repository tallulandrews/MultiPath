#' Preprocessed Gene Ontology: GOBP_layers
#'
#' @description Reorganization of the GO:BP ontology into a stack of layers at different depths from the root term. Individual GO terms appear in multiple layers if there are multiple routes of different length from the root to that term. Created: May 2025
#'
#' @format a list of 18 vectors:
#' \describe{
#'    \item{item name}{layer within the GO:BP ontology starting with 1 = root.}
#'    \item{value}{a vector of GO term ids for the terms that are reachable in exactly that many steps from the root.}
#'  }
#'
#' @source {GO.db} Bioconductor database package.
"GOBP_layers"

#' GTEx Tissue-specific Expression: gtex_genes
#'
#' @description Average expression (counts-per-million) as measured by bulkRNAseq in each tissue from the GTEx database (v10, RNASeqQCv2.4.2). Created: Jun 2025
#'
#' @format a matrix containing 20801 rows and 30 columns:
#' \describe{
#'    \item{columns}{correspond to 30 different tissues in alphabetical order: Adipose, Adrenal, Bladder, Blood, Blood Vessel, Brain, Breast, Cervix, Colon, Esophagus, Fallopian Tube, Heart, Kidney, Liver, Lung, Muscle, Nerve, Ovary, Pancreas, Pituitary, Prostate, Salivary, Skin, Small Intestine, Spleen, Stomach, Testis, Thyroid, Uterus, Vagina.}
#'    \item{rows}{genes, contains a mix of gene symbols and Ensembl IDs(if no gene symbol was available)}
#'  }
#'
#' @source \url{https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_reads.gct.gz}
"gtex_genes"


### Load info ###

#' Load the GOBP ontology annotation
#'
#' This load the GO Biological Process annotations for all human genes using the org.hs.eg.db package.
#'
#' @return a table with the following columnd: the ENTREZ gene id, the go term annotation, the evidence code, the Ontology (BP), and the gene symbol (GENE). 
#' @examples
#' gene2gobp <- load_GOBP_info()
#' @export
load_GOBP_info <- function() {
	#library('org.Hs.eg.db')
	gene2gobp <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, columns=c("GO", "ONTOLOGY"), 
					keys=AnnotationDbi::keys(org.Hs.eg.db::org.Hs.eg.db))
	gene2gobp <- gene2gobp[gene2gobp$ONTOLOGY == "BP",]
	entrez2symbol <- org.Hs.eg.db::org.Hs.egSYMBOL
	mapped_genes <- AnnotationDbi::mappedkeys(entrez2symbol)
	entrez2symbol <- as.list(entrez2symbol[mapped_genes])
	gene2gobp <- gene2gobp[gene2gobp$ENTREZID %in% names(entrez2symbol),]
	gene2gobp$GENE <- unlist(entrez2symbol[gene2gobp$ENTREZID])
	return(gene2gobp)
}

#' Convert Gene2GO into pathways list
#'
#' Converts the GOBP annotations used for simulations into appropriate format for use in do_ora and do_gsea
#'
#' @param annotation_table table of gene annotations from 'load_GOBP_info'
#' @param gene_col the name of the column of the table with geneID
#' @param anno_col the name of the column of the table with the GO term ID
#'
#' @return a list of GO pathways with their associated genes
#' @export
convert_GOBP <- function(annotation_table, gene_col="GENE", anno_col="GO") {
	pathways <- list()
	for (term in unique(annotation_table[,anno_col])) {
		pathways[[term]] <- annotation_table[ annotation_table[,anno_col] == term ,gene_col]	
	}
	return(pathways)
}

### Get Tissue Specific Genes ###

#' Obtain tissue-specific genes
#'
#' Get the tissue specific genes for a random selection of tissues.
#'
#' Randomly selects "n" tissues and gets the a list of genes with an average expression > threshold counts per million (cpm).
#' 
#' @param n the number of tissues to randomly sample.
#' @param seed set a random seed for the sampling.
#' @param threshold the cpm threshold to define a gene as expressed in the tissue.
#'
#' @return a list of tissue specific genes for each selected tissue.
#' @examples
#' genes <- get_tissue_specific_genes(3, seed=3819, threshold=5)
#' @export
get_tissue_specific_genes <- function(n, seed=NULL, threshold=5){
	if (!is.null(seed)) {
		set.seed(seed)
	}
	#load("gtex_genes.rda")
	tissues <- sample(1:ncol(gtex_genes), size=n, replace=TRUE)
	genes <- lapply(tissues, function(t){
		rownames(gtex_genes)[gtex_genes[,t] > threshold]
	})
	names(genes) <- colnames(gtex_genes)[tissues]
	return(genes)
}

### Draw Genes ###

#' Draw Genes from Pathways
#'
#' Randomly draws true positive and false positive genes from a set of provided pathways.
#'
#' Given a set of pathways and a set of tissue-specific genes, this function randomly selects a set number of genes from each pathway that are also tissue-specific. Subsequently randomly selects genes uniformly from the tissue-specific genes as false-positive genes. Warning: no error checking / input validation is preformed for this function.
#' 
#' @param gt_pathways a vector of pathways to draw genes from.
#' @param N the number of genes to draw from each pathway.
#' @param fpr the proportion of the returned genes that should be false-positives.
#' @param genes a vector of gene IDs that constitute the background - i.e. all detectable genes.
#' @param gene_anno a table of gene-pathway annotations, must have a column named "GO" containing the pathway term, and a column named "GENE" with the gene ids.
#' @param seed the random seed.
#'
#' @return a vector of selected genes.
#' 
#' @examples
#' gene2go <- load_GOBP_info()
#' tissue_specific <- get_tissue_specific_genes(1, seed=301)
#' genes <- draw_genes(c("GO:1990573","GO:0051898"), N=5, fpr=0.5, tissue_specific[[1]], gene2go, seed=101)
#' genes <- draw_genes(c("GO:1990573","GO:0051898"), N=5, fpr=0, tissue_specific[[1]], gene2go, seed=101)
#' @export
draw_genes <- function(gt_pathways, N=50, fpr=0.05, genes, gene_anno, seed=NULL){
	if (!is.null(seed)) {
		set.seed(seed)
	}
	this_genes <- c();
	gene_anno <- gene_anno[gene_anno[,"GENE"] %in% unlist(genes),]
	for (term in gt_pathways){
		tmp_genes <- sample(gene_anno[gene_anno[,"GO"]==term,"GENE"], size=N)
		this_genes <- c(this_genes, tmp_genes)
	}
	this_genes <- unique(this_genes)
	# False Positive Genes
	n_FP <- ceiling(fpr*length(this_genes) / (1-fpr))
	unselected <- unique(gene_anno[,"GENE"])
	unselected <- unselected[unselected %in% unlist(genes)] 
	unselected <- unselected[!unselected %in% this_genes]
	FP_genes <- sample(unselected, size=n_FP)
	return(unique(c(this_genes, FP_genes)))
}

### Generate log2fc ###

#' Generate generate log2fcs
#'
#' Generates example log2fcs from GTEX data.
#'
#' Randomly selects two tissues from the GTEX database, and calculates the log2fc for every gene betweeen them. Excludes any NAs or infinite values resulting from a gene being 0 in the denominator.
#'
#' @param seed the random seed.
#' @returns a vector of log2 fold changes.
#' @examples
#' l2fc <- generate_log2fcs(seed=1938)
#' @export
generate_log2fcs <- function(seed=NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}
	# generate realistic log2fc from GTEX
	#load("gtex_genes.rda")
	tmp <- sample(1:ncol(gtex_genes),size=2)
	log2fcs <- log2(gtex_genes[,tmp[1]]/gtex_genes[,tmp[2]])
	log2fcs <- log2fcs[!is.na(log2fcs)]
	log2fcs <- log2fcs[is.finite(log2fcs)]
	return(log2fcs)
}


### Generate GSEA Input ###

#' Generate GSEA input
#'
#' Generates a vector of log2fcs for GSEA with a set of true genes at the top.
#'
#' Creates a vector of weighted genes using provided log2fcs or generating these values using 'generate_log2fcs' if not provided. Then assigns the top highest positive values to the provided 'true_genes' after randomly permuting them, the rest of the values are assigned to the background genes provided in 'genes' which are also permuted to ensure no possible patterns. If there are insufficient log2fcs provided, the missing values are filled in with 0s.
#'
#' @param true_genes true-positive genes to be placed at the top.
#' @param genes background genes.
#' @param log2fcs vector of log2fcs, if NULL log2fcs are generated using 'generate_log2fcs'
#' @param seed the random seed to set.
#' 
#' @return a named vector of log2fcs with the true genes at the top (most positive log2fcs).
#'
#' @examples
#' tissue_specific <- get_tissue_specific_genes(1, seed=301)
#' TP_genes <- sample(tissue_specific[[1]], 50)
#' gsea_input <- generate_gsea_input(TP_genes, tissue_specific[[1]], seed=2911) 
#' gsea_input <- generate_gsea_input(TP_genes, tissue_specific[[1]], log2fcs=rnorm(2000), seed=2911) 
#' @export
generate_gsea_input <- function(true_genes, genes, log2fcs=NULL, seed=NULL) {
	if (!is.null(seed)) {
		set.seed(seed)
	}

	genes <- unlist(genes)
	if (is.null(log2fcs)) {
		log2fcs <- generate_log2fcs()
	}
	log2fcs <- log2fcs[1:length(genes)]
	log2fcs[is.na(log2fcs)] <- 0
	log2fcs <- log2fcs[order(log2fcs, decreasing=TRUE)]
	true_genes <- sample(true_genes, size=length(true_genes))
	other_genes <- genes[!genes %in% true_genes]
	other_genes <- sample(other_genes, size=length(other_genes))
	new_names <- c(true_genes, other_genes)
	names(log2fcs) <- new_names
	return(log2fcs)
}

### Selecting Pathways Functions ###

#' Random Pathways
#' 
#' Select random tissue-specific pathways.
#'
#' Using various contraints randomly samples valid GO pathways.
#'
#' @param n the number of GO terms to select
#' @param N minimum number of genes per pathway.
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param seed the random seed to ensure reproducibility
#' @param genes a background set of genes to limit pathways to.
#' @param gene_anno table of gene->GOterm annotations (from: load_GOBP_info())
#'
#' @return a vector of GO term IDs fitting the requirements.
#' @examples
#' gene2go <- load_GOBP_info()
#' tissue_specific <- get_tissue_specific_genes(1, seed=301)
#' terms <- select_random_terms(5, N=20, max_term_size=1000, seed=2938, genes=tissue_specific, gene_anno=gene2go)
#' @import GO.db
#' @export
 
select_random_terms <- function(n, N=50, max_term_size=1000, seed=2938, genes=NULL, gene_anno) {
	set.seed(seed)
	#requireNamespace('GO.db')
	# Subset to viable pathways give the provided genes and N
	if (!is.null(genes)) {
		gene_anno <- gene_anno[gene_anno$GENE %in% unlist(genes),]
	}
	pathway_size <- table(gene_anno$GO)
	viable_pathways <- names(pathway_size)[pathway_size > N & pathway_size < max_term_size]
	return(sample(viable_pathways, n, replace=FALSE))
}

### Hierarchical Terms ###

#' Hierarchical Pathways
#' 
#' Select random pairs of pathways separated in the GO hierarchy by a specified number of steps.
#'
#' Uses pre-calculated lists of terms at each depth in the GOBP ontology (using GO.db package for ontology structure) to randomly select 'n' pathways. Then walks up the ontology the specified 'depth' number of steps to select the second pathway. Various additional constraints can be imposed on the selected pathways. This function makes 'n_attempts' attempts to find pairs of pathways that fit all the constraints.
#'
#' @param depth the number of steps through the ontology separating each pair of pathways.
#' @param n the number of GO term pairs to select
#' @param N minimum number of genes per pathway.
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param seed the random seed to ensure reproducibility
#' @param genes a background set of genes to limit pathways to.
#' @param gene_anno table of gene->GOterm annotations (from: load_GOBP_info())
#' @param n_attempts the maximum number of attempts to make to find pathways that fit all the criteria.
#'
#' @return a vector of GO term IDs, pairs are sequential: 1,2,1,2,1,2,1,2,...
#' @examples
#' gene2go <- load_GOBP_info()
#' tissue_specific <- get_tissue_specific_genes(1, seed=301)
#' terms <- select_hierarchy_term_pairs(depth=1, 2, N=10, max_term_size=1000, seed=2938, genes=tissue_specific, gene_anno=gene2go)
#' terms <- select_hierarchy_term_pairs(depth=5, 2, N=10, max_term_size=1000, seed=2938, genes=tissue_specific, gene_anno=gene2go)
#' child_terms = terms[seq(from=1, to=length(terms)-1, by=2)]
#' parent_terms = terms[seq(from=2, to=length(terms), by=2)]
#' @import GO.db
#' @export
select_hierarchy_term_pairs <- function(depth, n, N=50, max_term_size=1000, seed=2938, genes=NULL, gene_anno, n_attempts=100){
	set.seed(seed)
	#requireNamespace('GO.db')
	#load('GOBP_layers.rda')
	term2parents <- as.list(GOBPPARENTS)
	# Read in necessary info:
	
	# Subset to viable pathways give the provided genes and N
	if (!is.null(genes)) {
		gene_anno <- gene_anno[gene_anno$GENE %in% unlist(genes),]
	}
	pathway_size <- table(gene_anno$GO)
	viable_pathways <- names(pathway_size)[pathway_size > N & pathway_size < max_term_size]
	pathway_layers <- intersectToList(GOBP_layers, viable_pathways)
	
	chosen <- c();
	check <- 0
	while(length(unique(chosen)) < 2*n & check < n_attempts) {
		possible <- unlist(pathway_layers[depth:length(pathway_layers)])
		possible <- possible[!possible %in% chosen]
		pick <- sample(possible, 1)
		current <- pick
		previous <- pick
		for (d in 1:depth) {
			parents <- c()
			previous <- c(previous, current);
			for (t in current) {
				parents <- c(parents, term2parents[[t]])
			}
			if (identical(current, parents)| length(parents) == 0) {
				break;
			} else {
				current <- parents
			}
		}

		current <- current[!current %in% previous]
		current <- current[current %in% viable_pathways]

		if (d == depth & !identical(current, pick) & length(current) > 0){
			p <- sample(current, 1)
			chosen <- c(chosen, c(pick, p)) 
		}
		check <- check+1
	}

	if (check == n_attempts) {
		warning("Warning: Ran out of attempts to find pathways fitting the criteria. Returning only those that were found rather than the number requested.")
	}
	return(chosen)
}

### Shared Pathways ###

#' Select Shared Pathways
#'
#' Generates a set of shared and unique pathways for modelling a multiple comparison
#'
#' Simulates a situation where multiple comparisons have been performed for one experiment and pathway analysis has been performed for each comparison. It assumes some pathways are shared amongst all comparisons and some are unique to each individual comparison. Unique pathways are drawn independently for each comparison thus may overlap but are always distinct from the shared pathways.
#' 
#' @param nc the number of comparisons
#' @param ns the number of pathways shared amongst all comparisons
#' @param nu the number of pathways unique to a comparison
#' @param N minimum number of genes per pathway.
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param seed the random seed to ensure reproducibility
#' @param genes a background set of genes to limit pathways to.
#' @param gene_anno table of gene->GOterm annotations (from: load_GOBP_info())
#'
#' @return a list of tables of pathways and whether or not the pathway is shared or unique.
#'
#' @examples
#' gene2go <- load_GOBP_info()
#' tissue_specific <- get_tissue_specific_genes(1, seed=301)
#' term_list <- select_shared_pathways(3, 2, 5, N=20, max_term_size=1000, seed=2938, genes=tissue_specific, gene_anno=gene2go)
#' term_list[[1]] # pathways for comparison 1
#' shared <- term_list[[1]][ term_list[[1]][,2] == "shared" ,1]
#' unique1 <- term_list[[1]][ term_list[[1]][,2] == "unique" ,1]
#' @export

select_shared_pathways <- function(nc, ns, nu, N=50, max_term_size=1000, seed=2938, genes=NULL, gene_anno) {
	set.seed(seed)
	#library('GO.db')
	# Subset to viable pathways give the provided genes and N
	if (!is.null(genes)) {
		gene_anno <- gene_anno[gene_anno$GENE %in% unlist(genes),]
	}
	pathway_size <- table(gene_anno$GO)
	viable_pathways <- names(pathway_size)[pathway_size > N & pathway_size < max_term_size]

	# Shared pathways
	shared_paths <- sample(viable_pathways, ns, replace=FALSE)
	viable_pathways2 <- viable_pathways[!viable_pathways %in% shared_paths]
	output <- list()
	for( c in 1:nc) {
		output[[c]] <- cbind(c(shared_paths, sample(viable_pathways2, nu, replace=FALSE)), 
				     c(rep("shared", length(shared_paths)), rep("unique", nu)))
		viable_pathways2 <- viable_pathways2[!viable_pathways2 %in% output[[c]]]
	}
	return(output)
}


#' Generate Random Pathways Simulation
#'
#' Generates DE gene output for ORA or GSEA for randomly selected pathways.
#'
#' @details Combines multiple functions to first randomly select a tissue for each simulation ('get_tissue_specific_genes'), select GO:BP terms ('select_random_terms'), randomly draw genes from those terms ('draw_genes'), generate log2fcs which are identical for all replicates ('generate_log2fcs'), and create the gsea input ('generate_gsea_input'). 
#'
#'This function is for the random pathways case, where terms are drawn at random equally likely across all valid pathways. 
#'
#' Constraints imposed on pathways are: tissue-specific gene expression, sufficient genes for the number of genes per pathway required, and no more than 'max_term_size' genes. Gene number constraints are imposed after filtering for tissue-specific gene expression.
#' 
#' @param n_replicas number of simulations to generate
#' @param n_path number of pathways selected per simulation
#' @param n_genes_per_path number of genes significantly DE for each selected pathway
#' @param fpr proportion of DE genes that are false-positives 
#' @param seed the random seed to ensure reproducibility
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param gene2terms table of gene->GOterm annotations (from: load_GOBP_info())
#'
#' @return a list of simulation output, each simulation output is a list that contains:
#'	gt_terms : the truly DE GO terms (vector)
#'	de_genes : the DE genes including false-positives (vector), needed for do_ora
#'	background : the tissue-specific genes used as the background (vector), needed for do_ora
#'	gsea_input : named vector of fold changes with the de_genes at the top, needed for do_gsea
#'
#' @examples
#' gene2go <- load_GOBP_info()[1:50000,] # reduced to limit runtime
#' sims <- create_rand_pathway_sim(1, 5, 15, fpr=0.05, gene2terms=gene2go)
#' paths <- convert_GOBP(gene2go)
#' ora_out <- do_ora(sims[[1]]$de_genes,pathways=paths, background=sims[[1]]$background)
#' gsea_out <- do_gsea(sims[[1]]$gsea_input, pathways=paths, nperm=100) # reduced number of permutations to limit runtime 
#' @export
create_rand_pathway_sim <- function(n_replicas, n_path, n_genes_per_path, fpr, seed=2938, max_term_size=1000, gene2terms) {	
	set.seed(seed)
	future_seed <- ceiling(stats::runif(100)*10000) + 1000
	seed_i = 1;

	# Pick a tissue for each replica
	tissue_genes <- get_tissue_specific_genes(n_replicas, seed=2903, threshold=5)
	
	log2fcs <- generate_log2fcs()
	out <- list()
	this = 1;
	for( genes in tissue_genes ) {
		this_name <- paste("rand",n_path, n_genes_per_path, fpr, this, sep="_");
		this <- this+1
		terms <- select_random_terms(n=n_path, N=n_genes_per_path, max_term_size=max_term_size, 
				gene_anno=gene2terms, genes=genes, seed=future_seed[seed_i]); seed_i <- seed_i+1;
		
		de_genes <- draw_genes(gt_pathways=terms, N=n_genes_per_path, fpr=fpr, genes=genes, 
				gene_anno=gene2terms, seed=future_seed[seed_i]); seed_i <- seed_i+1;
		
		fgsea_input <- generate_gsea_input(true_genes=de_genes, genes, log2fcs=log2fcs)
		sim <- list(gt_terms = terms,  # ground truth
				de_genes = de_genes, background=unlist(genes), # input for ORA
				gsea_input = fgsea_input) # input of gsea
		out[[this_name]] <- sim
	}
	return(out)
}

#' Generate Hierarchical Pathways Simulation
#'
#' Generates DE gene output for ORA or GSEA for pairs of related pathways.
#'
#' @details Combines multiple functions to first randomly select a tissue for each simulation ('get_tissue_specific_genes'), select GO:BP terms ('select_hierarchy_term_pairs'), randomly draw genes from those terms ('draw_genes'), generate log2fcs which are identical for all replicates ('generate_log2fcs'), and create the gsea input ('generate_gsea_input'). 
#'
#'This function is for the hierarchical pathways case, where pairs of terms are drawn such that they are separated from each other by 'depth' steps in the GOBP ontology structure. 
#'
#' Constraints imposed on pathways are: tissue-specific gene expression, sufficient genes for the number of genes per pathway required, and no more than 'max_term_size' genes. Gene number constraints are imposed after filtering for tissue-specific gene expression.
#' 
#' @param n_replicas number of simulations to generate
#' @param n_path number of pathways selected per simulation (should be an even number as pathways are selected as pairs)
#' @param n_genes_per_path number of genes significantly DE for each selected pathway
#' @param fpr proportion of DE genes that are false-positives 
#' @param depth the number of steps through the ontology separating each pair of pathways.
#' @param seed the random seed to ensure reproducibility
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param gene2terms table of gene->GOterm annotations (from: load_GOBP_info())
#'
#' @return a list of simulation output, each simulation output is a list that contains:
#'	gt_terms : the truly DE GO terms (vector)
#'	de_genes : the DE genes including false-positives (vector), needed for do_ora
#'	background : the tissue-specific genes used as the background (vector), needed for do_ora
#'	gsea_input : named vector of fold changes with the de_genes at the top, needed for do_gsea
#'
#' @examples
#' gene2go <- load_GOBP_info()[1:50000,] # reduced to limit runtime.
#' sims <- create_hierarchical_pathway_sim(1, 6, 15, depth=2, fpr=0.05, gene2terms=gene2go)
#' paths <- convert_GOBP(gene2go)
#' ora_out <- do_ora(sims[[1]]$de_genes, pathways=paths, background=sims[[1]]$background)
#' gsea_out <- do_gsea(sims[[1]]$gsea_input, pathways=paths, nperm=100) # reduced number of permutations to limit runtime 
#' @export
create_hierarchical_pathway_sim <- function(n_replicas, n_path, n_genes_per_path, fpr, depth, seed=2938, max_term_size=1000, gene2terms) {	
	set.seed(seed)
	future_seed <- ceiling(stats::runif(100)*10000) + 1000
	seed_i = 1;

	# Pick a tissue for each replica
	tissue_genes <- get_tissue_specific_genes(n_replicas, seed=2903, threshold=5)
	
	log2fcs <- generate_log2fcs()
	out <- list()
	this = 1;
	for( genes in tissue_genes ) {
		this_name <- paste("hier", depth, n_path, n_genes_per_path, fpr,this, sep="_");
		this <- this+1
		terms <- select_hierarchy_term_pairs(n=ceiling(n_path/2), N=n_genes_per_path, max_term_size=max_term_size, depth=depth, 
				gene_anno=gene2terms, genes=genes, seed=future_seed[seed_i]); seed_i <- seed_i+1;
		
		de_genes <- draw_genes(gt_pathways=terms, N=n_genes_per_path, fpr=fpr, genes=genes, 
				gene_anno=gene2terms, seed=future_seed[seed_i]); seed_i <- seed_i+1;
		
		fgsea_input <- generate_gsea_input(true_genes=de_genes, genes, log2fcs=log2fcs)
		sim <- list(gt_terms = terms,  # ground truth
				de_genes = de_genes, background=unlist(genes), # input for ORA
				gsea_input = fgsea_input) # input of gsea
		out[[this_name]] <- sim
	}
	return(out)
}


#' Generate Shared Pathways Simulation
#'
#' Generates DE gene output for ORA or GSEA for a multiple comparison experiment.
#'
#' @details Combines multiple functions to first randomly select a tissue for each simulation ('get_tissue_specific_genes'), select GO:BP terms ('select_shared_pathways'), randomly draw genes from those terms ('draw_genes'), generate log2fcs which are identical for all replicates ('generate_log2fcs'), and create the gsea input ('generate_gsea_input'). 
#'
#'This function is for the shared case, where multiple sets of terms/genes are drawn representing different comparisons of the same experiment. Some pathways are shared amongst all comparisons, some are not - i.e. unique to a comparison. 
#'
#' Constraints imposed on pathways are: tissue-specific gene expression, sufficient genes for the number of genes per pathway required, and no more than 'max_term_size' genes. Gene number constraints are imposed after filtering for tissue-specific gene expression.
#' 
#' @param n_replicas number of simulations to generate
#' @param n_comparisons number of comparisons in the simulation 
#' @param n_path_shared number of pathways shared amongst all comparison in each simulation
#' @param n_path_unique number of pathways unique to each comparison in each simulation
#' @param n_genes_per_path number of genes significantly DE for each selected pathway
#' @param fpr proportion of DE genes that are false-positives 
#' @param seed the random seed to ensure reproducibility
#' @param max_term_size the maximum genes per pathway to consider that pathway
#' @param gene2terms table of gene->GOterm annotations (from: load_GOBP_info())
#'
#' @return a list of simulation output, each simulation output is a list that contains:
#'	gt_terms : the truly DE GO terms (list for each comparison - includes a column indicating if the pathway is shared or unique)
#'	de_genes : the DE genes including false-positives (list for each comparison), needed for do_ora
#'	background : the tissue-specific genes used as the background (vector), needed for do_ora
#'	gsea_input : named vector of fold changes with the de_genes at the top (list for each comparison), needed for do_gsea
#'
#' @examples
#' gene2go <- load_GOBP_info()
#' sims <- create_rand_pathway_sim(1, 5, 15, fpr=0.05, gene2terms=gene2go)
#' paths <- convert_GOBP(gene2go)
#' ora_out <- do_ora(sims[[1]]$de_genes, pathways=paths, background=sims[[1]]$background)
#' gsea_out <- do_gsea(sims[[1]]$gsea_input, pathways=paths, nperm=100) # reduced number of permutations to limit runtime 
#' @export
create_shared_pathway_sim <- function(n_replicas, n_comparisons = 5, n_path_shared, n_path_unique,  n_genes_per_path, fpr,  seed=2938, max_term_size=1000, gene2terms) {	
	set.seed(seed)
	future_seed <- ceiling(stats::runif(100)*10000) + 1000
	seed_i = 1;

	# Pick a tissue for each replica
	tissue_genes <- get_tissue_specific_genes(n_replicas, seed=2903, threshold=5)
	
	log2fcs <- generate_log2fcs()
	out <- list()
	this = 1;
	for( genes in tissue_genes ) {
		this_name <- paste("shared", n_path_shared, n_path_unique, n_genes_per_path, fpr,this, sep="_");
		this <- this+1
		terms <- select_shared_pathways(nc=n_comparisons, ns=n_path_shared, nu=n_path_unique, 
				N=n_genes_per_path, max_term_size=max_term_size,
				gene_anno=gene2terms, genes=genes, seed=future_seed[seed_i]); seed_i <- seed_i+1;

		de_genes <- lapply(terms, function(x) {draw_genes(x[,1], N=n_genes_per_path, fpr=fpr, genes=genes,
                                gene_anno=gene2terms, seed=NULL)})
		names(de_genes) <- paste0("comparison", 1:length(de_genes))

		fgsea_input <- lapply(de_genes, generate_gsea_input, genes=genes, log2fcs=log2fcs)		

		sim <- list(gt_terms = terms,  # ground truth - I need to separate this into shared & unique gt.
				de_genes = de_genes, background=unlist(genes), # input for ORA
				gsea_input = fgsea_input) # input of gsea
		out[[this_name]] <- sim
	}
	return(out)
}
