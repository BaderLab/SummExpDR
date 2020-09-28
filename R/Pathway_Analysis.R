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

#' Fisher's Exact Test on all pathways in gmt given query geneset
#' @param query character vector with query feature set
#' @param gmt_file gmt file containing pathways of interest
#' @param FDR FDR to filter at
#' @param n_cores
#' @param pbar
#' @value data frame denoting pathway name, overlap, overlapping genes, and number of genes in pathway present in master set. Also includes p-value + FDR for overlap.
#' @export

do_fisher_multi <- function(query, gmt_file, master_set, FDR = 0.10, n_cores = 1, pbar = FALSE) {
  capt_out <- capture.output(GSA_gmt_output <- GSA::GSA.read.gmt(gmt_file))
  pathway_list <- GSA_gmt_output$genesets
  names(pathway_list) <- GSA_gmt_output$geneset.names
  do_fisher_single <- function(x) {
    pathway_x <- pathway_list[[x]]
    pathway_x_subs <- intersect(pathway_x, master_set)
    fisher_test_result <- fisher.test.sets(set1 = query, set2 = pathway_x_subs, master_set = master_set)
    intersection <- intersect(query, pathway_x_subs)
    return(data.frame(pathway_name = x,
                      intersection = paste(intersection, collapse = ';'),
                      pathway_kept = paste(pathway_x_subs, collapse = ';'),
                      pathway_orig = paste(pathway_x, collapse = ';'),
                      n_intersect = length(intersection),
                      n_kept = length(pathway_x_subs),
                      n_orig = length(pathway_x),
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
