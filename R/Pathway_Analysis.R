


#' Load gmt file
#' @param gmt_file gmt_file to load in
#' @param min_genes minimum genes per pathway. Note that min filter is applied after filtering for genes in master set
#' @param max_genes maximum genes per pathway. Note that max filter is applied after filtering for genes in master set
#' @param master_set master gene set, of which all pathways must be a subset in the result
#' @param n_cores number of cores to use for filtering
#' @param pbar whether or not to use a progressbar
#' @details in addition to filtering for genes (strings) in master set, filters out empty characters (\'\').
#' @value named list of pathways (character vectors) filtered according to criteria set in parameters
#' @export

load_gmt <- function(gmt_file, min_genes = 15, max_genes = 200, master_set = NULL, n_cores = 1, pbar = FALSE) {
  capt_out <- capture.output(GSA_gmt_output <- GSA::GSA.read.gmt(gmt_file))
  pathway_list <- GSA_gmt_output$genesets
  pathway_names <- GSA_gmt_output$geneset.names
  if (is.null(master_set)) {
    master_set <- unique(unlist(pathway_list))
  }
  character_filter <- function(x, master_set) {
    in_master <- x %in% master_set
    if (sum(in_master) < 1) {
      x <- character(0)
    } else {
      x <- x[in_master]
    }
    not_empty <- x != ''
    if (sum(not_empty) > 0) {
      x <- x[not_empty]
    } else {
      x <- character(0)
    }
    return(x)
  }
  size_filter <- function(x, min_genes, max_genes) {
    len_x <- length(x)
    return((len_x >= min_genes) && (len_x <= max_genes))
  }
  if (n_cores > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  }
  num_pathways <- length(pathway_list)
  print(paste('Filtering genes in', num_pathways, 'pathways'))
  pathway_list <- BiocParallel::bplapply(pathway_list, FUN = function(x) {character_filter(x, master_set)}, BPPARAM = bpparam)
  names(pathway_list) <- pathway_names
  print(paste('filtering pathways for specified size'))
  keep_pathways <- unlist(BiocParallel::bplapply(pathway_list, FUN = function(x) {size_filter(x, min_genes, max_genes)}, BPPARAM = bpparam))
  num_keep <- sum(keep_pathways)
  print(paste('kept', num_keep, 'pathways'))
  pathway_list <- pathway_list[keep_pathways]
  return(pathway_list)
}

#' Fisher's Exact Test on all pathways in gmt given query geneset
#' @param query character vector with query feature set
#' @param pathway_list list of pathways of interest
#' @param FDR FDR to filter at
#' @param n_cores number of cores to use for fisher's exact test
#' @param pbar whether or not to use a progressbar
#' @value data frame denoting pathway name, overlap, overlapping genes, and number of genes in pathway present in master set. Also includes p-value + FDR for overlap.
#' @export

do_fisher_multi <- function(query, pathway_list, master_set, FDR = 0.10, n_cores = 1, pbar = FALSE) {
  stopifnot(!is.null(names(pathway_list)))
  do_fisher_single <- function(x) {
    pathway_x <- pathway_list[[x]]
    fisher_test_result <- fisher_test_sets(set1 = query, set2 = pathway_x, master_set = master_set)
    intersection <- intersect(query, pathway_x)
    return(data.frame(pathway_name = x,
                      intersection = paste(intersection, collapse = ';'),
                      pathway = paste(pathway_x, collapse = ';'),
                      n_intersect = length(intersection),
                      n_pathway = length(pathway_x),
                      p_val = fisher_test_result$p.value,
                      stringsAsFactors = FALSE)
           )
  }
  if (n_cores == 1) {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  }
  df_list <- BiocParallel::bplapply(names(pathway_list), BPPARAM = bpparam, FUN = function(x) {do_fisher_single(x)})
  result_df <- as.data.frame(dplyr::bind_rows(df_list))
  result_df$FDR <- p.adjust(result_df$p_val, method = 'BH')
  result_df <- result_df[result_df$FDR <= FDR, ]
  return(result_df)
}

#' Run GseaPreranked from R via system call
#'
#' @param gsea_rank_list_path  character, length 1. path to ranked list file
#' @param gsea_analysis_name character, length 1. name for gsea analysis
#' @param gsea_rpt_label character, length 1. report label
#' @param gsea_jar_path path to GSEA jar file
#' @param gsea_out directory under which GSEA results will be output. GSEA makes a directory
#' under this one to contain the GSEA output
#' @param gsea_gmx path to gmt file
#' @param gsea_num_permutations number of permutations to run for GSEA
#' @param gsea_min_gs_size minimum geneset size for GSEA
#' @param gsea_max_gs_size maximum geneset size for GSEA
#' @param plot_top generate plots for top (plot_top) positively and negatively enriched pathways
#' @param seed numeric, length 1. seed to use
#' @param mem character, length 1. memory (in GB) to allot to GSEA. E.g. 16G for 16 GB
#' @details assumes you have GSEA version 3.0 .jar file in your filesystem
#'
#' @return exit status of system call to GseaPreranked , numeric of length 1.
#' @export
#'
#' @examples
enrich_GSEA <- function(gsea_rank_list_path, # path to ranked list file
                        gsea_analysis_name = 'analysis',
                        gsea_rpt_label = 'report',
                        gsea_jar_path = file.path("~/GSEA/gsea-3.0.jar"),
                        gsea_out = NULL, # specify as null or as a string denoting valid filepath
                        gsea_gmx = NULL, # specify as null or as a string denoting path to .gmt file
                        gsea_num_permutations = 1000,
                        gsea_min_gs_size = 15,
                        gsea_max_gs_size = 200,
                        plot_top = 20,
                        seed = 12345,
                        mem = '16G' # max memory to allocate to Java program
) {
  print('enrich_GSEA')
  print(Sys.time())
  ### Declare user-defined settings
  stopifnot(is.character(gsea_jar_path))
  if (!file.exists(gsea_jar_path)) {
    stop(paste('cannot find file', gsea_jar_path))
  }
  if (is.null(gsea_out)) {
    gsea_out <- file.path(getwd(), "gsea_output")
  }

  ## check for existing directory
  if (!dir.exists(gsea_out)) {
    stop(paste('output directory', gsea_out, 'does not exist'))
  }

  if(is.null(gsea_gmx)) {

    ## Get most up to date human GMT file from baderlab.org genesets database.
    ## Below code for getting current baderlab gmt file essentially ripped off
    ## of Ruth Isserlin's code from Enrichment Map pipeline:
    ## https://baderlab.github.io/Cytoscape_workflows/EnrichmentMapPipeline/Protocol2_createEM.html
    ## code was copied Dec 1, 2018
    gmt_url = "http://download.baderlab.org/EM_Genesets/current_release/Human/symbol/"

    #list all the files on the server
    filenames = RCurl::getURL(gmt_url)
    tc = textConnection(filenames)
    contents = readLines(tc)
    close(tc)

    #get the gmt that has all the pathways and does not include terms inferred from electronic annotations(IEA)
    #start with gmt file that has pathways only
    rx = gregexpr("(?<=<a href=\")(.*.GOBP_AllPathways_no_GO_iea.*.)(.gmt)(?=\">)",
                  contents, perl = TRUE)
    gmt_file = unlist(regmatches(contents, rx))

    dest_gmt_file <- file.path(gsea_out, gmt_file, sep="")
    dest_gmt_file <- sub('.gmt/$', '.gmt', dest_gmt_file)
    download.file(
      paste(gmt_url,gmt_file,sep=""),
      destfile=dest_gmt_file
    )
    gsea_gmx <- dest_gmt_file
  } else if (is.character(gsea_gmx)) {
    if (!file.exists(gsea_gmx)) {
      stop(paste('could not find file', gsea_gmx))
    }
  } else {
    stop(paste('gsea_gmx must be a string denoting valid filepath'))
  }


  ## make command
  command <- paste("java -cp", gsea_jar_path,
                   paste0("-Xmx",mem," xtools.gsea.GseaPreranked"),
                   "-rpt_label", gsea_analysis_name,
                   "-out", gsea_out,
                   "-gmx", gsea_gmx,
                   "-rnk", gsea_rank_list_path,
                   "-nperm", gsea_num_permutations,
                   "-set_min", gsea_min_gs_size,
                   "-set_max", gsea_max_gs_size,
                   "-collapse false",
                   "-scoring_scheme weighted",
                   "-permute gene_set",
                   "-num 100",
                   paste("-plot_top_x", as.character(plot_top)),
                   paste("-rnd_seed", as.character(seed)),
                   "-zip_report false",
                   "-gui false",
                   ">",
                   file.path(gsea_out, paste(gsea_analysis_name,'_', gsea_rpt_label, ".txt", sep="")),
                   sep=" ")
  ## EXECUTE COMMAND
  print(paste('using gmt file:', gsea_gmx))
  print(paste('executing command:', command))
  exit_status <- system(command)
  if (exit_status == 0) {
    print('enrich_gsea system call finished.')
  } else {
    print(paste('enrich_GSEA system call finished with exit status', exit_status))
    warning('enrich_GSEA system call finished with non-zero exit status')
  }

  print(Sys.time())
  return(exit_status)
}

### Read GSEA results
#' Extract results tables for positively and negatively enriched pathways
#'
#' @param gsea_results_dir character, length 1. path to GSEA results.
#'
#' @return list with 2 entries: na_pos_report = contents of gsea_report_for_na_pos\[0-9\]*.xls,
#' na_neg_report = contents of gsea_report_for_na_neg\[0-9\]*.xls
#' @export
#'
#' @examples
retrieve_gsea <- function(gsea_results_dir) {
  print('retrieve_gsea')
  pos_pathway_file <- dir(gsea_results_dir)[grep('gsea_report_for_na_pos_[0-9]*.xls',dir(gsea_results_dir))]
  neg_pathway_file <- dir(gsea_results_dir)[grep('gsea_report_for_na_neg_[0-9]*.xls',dir(gsea_results_dir))]
  if (length(pos_pathway_file) != 1) {
    print('pos pathway files')
    print(pos_pathway_file)
    stop(paste('found', as.character(length(pos_pathway_file)), 'matches for gsea_report_for_na_pos_[0-9]*.xls'))
  } else if (length(neg_pathway_file) != 1) {
    print('neg pathway files')
    print(neg_pathway_file)
    stop(paste('found', as.character(length(neg_pathway_file)), 'matches for gsea_report_for_na_neg_[0-9]*.xls'))
  }
  # load in the geneset files
  print(paste('retrieving in', gsea_results_dir))
  print(pos_pathway_file)
  print(neg_pathway_file)
  na_pos_report <- read.delim(file = file.path(gsea_results_dir, pos_pathway_file)) # positively enriched genesets
  na_neg_report <- read.delim(file = file.path(gsea_results_dir, neg_pathway_file)) # negatively enriched genesets
  GSEA_result <- list(na_pos_report = na_pos_report, na_neg_report = na_neg_report)
  return(GSEA_result)
}
