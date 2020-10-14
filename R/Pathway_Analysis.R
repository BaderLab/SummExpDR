#' Fisher's Exact Test, with set input
#'
#' @param set1 character vector containing a set of labels
#' @param set2 character vector containing a set of labels
#' @param master_set set of which set1 and set2 are a subset of and from which they were sampled
#' @description computes fisher's exact test for association between set1 and set2
#' given a master set from which set1 and set 2 were drawn.
#' @return the result of fisher.test run for a contingency table made given the input sets
#' @export
#'
#' @examples
fisher.test.sets <- function(set1, set2, master_set) {
  if (!all(set1 %in% master_set) || !all(set2 %in% master_set)) {
    stop('all members of set1 and set2 must be members of master set')
  }
  TT <- length(intersect(set1, set2))
  TF <- length(set2) - TT
  FT <- length(set1) - TT
  FF <- length(setdiff(master_set, set1)) - TF
  cont.mat <- matrix(c(TT, TF, FT, FF), nrow = 2)
  test.result <- fisher.test(cont.mat, alternative = 'greater')
  return(test.result)
}


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
    fisher_test_result <- fisher.test.sets(set1 = query, set2 = pathway_x, master_set = master_set)
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
