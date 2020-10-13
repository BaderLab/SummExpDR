###############################################################################
### Utilities for Analyzing DNAm Data ###

## Note that Illumina 850K (EPIC) array uses Gencode v12 annotations
##  https://www.gencodegenes.org/human/release_12.html
##  ../../data/raw/misc/infinium-methylationepic-manifest-column-headings.pdf
##  https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf

#' Map gene to probes
#'
#' @param gene = gene name
#' @param gene_probe_mapping = gene to probe mapping. must be dataframe containing genes in rows and column 'mapped_probes' of semicolon delimited strings
#' denoting EPIC probes
#' @value character vector of illumina EPIC array probes
#' @export
#'
get_mapped_probes <- function(gene, gene_probe_mapping) {
  if (!gene %in% rownames(gene_probe_mapping)) {
    stop(paste(gene, 'not found in mapping'))
  } else {
    probes <- strsplit(gene_probe_mapping[gene,'mapped_probes'], split = ';')[[1]]
  }
  return(probes)
}

#' Get check if given gene has min number of probes
#' @param min_probes = minimum number of probes for a gene
#' @param data_mat = nxm feature x sample matrix
#' @inheritParams get_mapped_probes
#' @value logical vector indicating which genes have minimum number of probes per gene
#' @export
check_min_probes <- function (min_probes, data_mat, gene_probe_mapping) {
  res <- as.logical(lapply(strsplit(gene_probe_mapping$mapped_probes, split = ';'),
                           FUN = function(x) {sum(x %in% rownames(data_mat)) >= min_probes}))
  names(res) <- gene_probe_mapping$gene
  if (!any(res)) {
    warning(paste('no gene has >=', min_probes, 'probes (min probes)'))
  }
  return(res)
}

#' Calculate Genewise Average Of Probe Feature
#'
#' @param inp_data = input data matrix, rows corresponding to probes, columns corresponding to features
#' @inheritParams get_mapped_probes
#' @param n_cores = number of cores to use
#' @return list with matrix of gene-wise averages, genes in rows, features in columns, and a second entry with number of probes per gene.
#' Note that genes without any probes are excluded from this matrix and probe count entry
#' @export
calculate_probe_avg <- function(inp_data, gene_probe_mapping, n_cores = 1, pbar = TRUE) {
  if (is.data.frame(inp_data)) {
    inp_data <- as.matrix(inp_data)
  } else if (!(is.matrix(inp_data) || is(inp_data, 'Matrix'))) {
    stop('must input a matrix or a data.frame (dataframe coerced to matrix)')
  }
  inp_data <- Matrix::Matrix(inp_data)
  probe_names <- rownames(inp_data)
  mapping <- strsplit(gene_probe_mapping$mapped_probes, split = ';')
  # keep track of genes whose probe mapping satisfies regex, has at least 1 probe
  probe_map_regex <- '^(([A-z0-9]|[.])+)+(;([A-z0-9]|[.])+)*$'
  has_unmatched_regex <- !grepl(probe_map_regex, gene_probe_mapping$gene)

  # set names of mapping to gene name prior to filtering
  names(mapping) <- gene_probe_mapping$gene
  # filter mapping
  genes_keep <- !has_unmatched_regex
  mapping <- mapping[genes_keep]
  n_gene <- length(mapping)
  do_avg <- function(i, map_ = mapping, mat_ = inp_data) {
    # return whether or not there was a successful mapping
    gene_name <- names(map_)[i]
    mapped_probes <- unlist(mapping[gene_name])
    probe_inds_i <- which(probe_names %in% mapped_probes)
    num_probes <- length(probe_inds_i)
    if (num_probes > 0) {
      data_subs_i = mat_[probe_inds_i, ]
      if (is.vector(data_subs_i)) {
        means_i <- data_subs_i
      } else {
        means_i = Matrix::colMeans(data_subs_i)
      }
    } else {
      # zero probes
      means_i <- rep(NA, ncol(mat_))
    }
    # give number of probes at end of vector
    means_i <- c(means_i, num_probes)
    # conversion to dataframe so that dplyr::bind_rows can handle
    means_i <- matrix(means_i, nrow = 1)
    colnames(means_i) <- 1:ncol(means_i)
    means_i <- as.data.frame(means_i)
    return(means_i)
  }

  if (n_cores > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  }
  loop_list <- BiocParallel::bplapply(1:n_gene, BPPARAM = bpparam, FUN = function(x) {do_avg(x)})
  loop_mat <- as.matrix(dplyr::bind_rows(loop_list))
  # extract probe avg mat, probe counts (probe counts are in last column)
  rownames(loop_mat) <- names(mapping)
  probe_avg_mat <- loop_mat[,1:(ncol(loop_mat) - 1)]
  colnames(probe_avg_mat) <- colnames(inp_data)
  probe_count <- loop_mat[,ncol(loop_mat)]
  names(probe_count) <- rownames(loop_mat)
  # count how many times you have unmatched regex or zero probes from a gene mapped to data's probes
  unmatched_regex_ct <- sum(has_unmatched_regex)
  has_zero_probes <- probe_count == 0
  zero_probes_ct <- sum(has_zero_probes)
  msg <- paste0(sum(probe_count), ' of ', length(probe_names), ' probes from dataset kept\n')
  msg <- paste0(msg, sum(genes_keep), ' of ', length(genes_keep), ' genes in gene-probe mapping kept\n')
  cat(msg)

  warn_msg <- ''

  if (unmatched_regex_ct > 0) {
    warn_msg <- paste0(warn_msg, unmatched_regex_ct, ' genes in provided mapping did not have proper regex for probe string\n')
  }
  if (zero_probes_ct > 0) {
    warn_msg <- paste0(warn_msg, zero_probes_ct, ' genes in provided mapping had 0 mapped probes in data\n')
  }
  if (warn_msg != '') {
    warning(warn_msg)
  }
  # filter for final probe_avg_mat and probe count vector to be returned by function
  probe_avg_mat <- probe_avg_mat[!has_zero_probes, ]
  probe_count <- probe_count[!has_zero_probes]

  return(list(probe_avg_mat = probe_avg_mat, probes_per_gene = probe_count))

}
