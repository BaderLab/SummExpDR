###############################################################################
### Utilities for Analyzing DNAm Data ###

## Note that Illumina 850K (EPIC) array uses Gencode v12 annotations
##  https://www.gencodegenes.org/human/release_12.html
##  ../../data/raw/misc/infinium-methylationepic-manifest-column-headings.pdf
##  https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/infinium-methylationepic-manifest-column-headings.pdf

#' Get probes mapped to by gene in data illumina_EPIC_hg19_10b4_gene_2_probe
#'
#' @param gene = gene name
#' @value character vector of illumina EPIC array probes
#' @export
#'
get_mapped_probes <- function(gene) {
  if (!gene %in% rownames(illumina_EPIC_hg19_10b4_gene_2_probe)) {
    stop(paste(gene, 'not found in '))
  } else {
    probes <- strsplit(illumina_EPIC_hg19_10b4_gene_2_probe[gene,'gene2probe'], split = ';')[[1]]
  }
  return(probes)
}

#' Get check if given gene has min number of probes
#' @param min_probes = minimum number of probes for a gene
#' @param data_mat = nxm feature x sample matrix
#' @value logical vector indicating which genes have minimum number of probes per gene
check_min_probes <- function (min_probes, data_mat) {
  res <- as.logical(lapply(strsplit(illumina_EPIC_hg19_10b4_gene_2_probe$gene2probe, split = ';'),
                           FUN = function(x) {sum(x %in% rownames(data_mat)) >= min_probes}))
  names(res) <- illumina_EPIC_hg19_10b4_gene_2_probe$gene_id
  if (!any(res)) {
    warning(paste('no gene has >=', min_probes, 'probes (min probes)'))
  }
  return(res)
}

#' Calculate Genewise Average Of Probe Feature
#'
#' @param inp_data = input data matrix, rows corresponding to probes, columns corresponding to features
#'
#' @return list with matrix of gene-wise averages, genes in rows, features in columns. Note that genes without any probes are excluded
#' from this matrix.
#'
#'
#' @export
calculate_probe_avg <- function(inp_data, mapping) {
  if (is.data.frame(inp_data)) {
    inp_data <- as.matrix(inp_data)
  } else if (!(is.matrix(inp_data) || is(inp_data, 'Matrix'))) {
    stop('must input a matrix or a data.frame (dataframe coerced to matrix)')
  }
  inp_data <- Matrix::Matrix(inp_data)
  # mapping <- illumina_EPIC_hg19_10b4_gene_2_probe$gene2probe
  # names(mapping) <- illumina_EPIC_hg19_10b4_gene_2_probe$gene_id
  # setup matrices/indexes
  n_gene <- length(mapping)
  probe_avg_mat <- matrix(data = numeric(n_gene*ncol(inp_data)), nrow = n_gene)
  rownames(probe_avg_mat) <- names(mapping)
  colnames(probe_avg_mat) <- colnames(inp_data)
  probes_per_gene <- integer(n_gene)
  names(probes_per_gene)
  probe_names <- rownames(inp_data)

  # keep track of genes whose probe mapping satisfies regex, has at least 1 probe
  probe_map_regex <- '^(([A-z0-9]|[.])+)+(;([A-z0-9]|[.])+)*$'
  unmatched_regex_ct <- 0
  zero_probes_ct <- 0
  genes_keep <- rep(TRUE, n_gene)
  probes_kept <- character(0)

  for (i in 1:n_gene) {
    gene_name <- names(mapping)[i]
    mapped_probe_str <- mapping[gene_name]
    if (!grepl(probe_map_regex, mapped_probe_str)) {
      unmatched_regex_ct <- unmatched_regex_ct + 1
      genes_keep[i] = FALSE
      next
    } else {
      probes_i <- strsplit(mapped_probe_str, split = ';')[[1]]
    }
    probe_inds_i <- probe_names %in% probes_i
    num_probes <- sum(probe_inds_i)
    if (num_probes > 0) {
      data_subs_i = inp_data[probe_inds_i, ]
      if (is.vector(data_subs_i)) {
        means_i <- data_subs_i
      } else {
        means_i = Matrix::colMeans(data_subs_i)
      }
      probe_avg_mat[i,] = means_i
      probes_kept <- union(probes_kept, probe_names[probe_inds_i])
      probes_per_gene[i] <- num_probes
    } else {
      genes_keep[i] <- FALSE
      zero_probes_ct <- zero_probes_ct + 1
      next
    }
  }

  msg <- paste0(length(probes_kept), ' of ', length(probe_names), ' probes from dataset kept\n')
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
  probe_avg_mat <- probe_avg_mat[genes_keep, ]
  probes_per_gene <- probes_per_gene[genes_keep]

  return(list(probe_avg_mat = probe_avg_mat, probes_per_gene = probes_per_gene))

}
