###############################################################################

#' Calculate correlation between all columns of matrices X and Y
#' @param X matrix
#' @param Y matrix
#' @param fdr_filter FDR filter for correlation values
#' @param self_only only calculate variables with same name in X and Y. example use would be
#' correlated DNA methylation for a gene and expression of same gene.
#' @value data.frame
#' @export

get_cor_vals <- function(X, Y, fdr_filter = 0.10, n_cores = 1, self_only = FALSE) {
  stopifnot(rownames(X) == rownames(Y))
  add_df <- function(x, y, loop_dataframe, p) {
    failed <- TRUE
    try({
      k <- loop_dataframe$i[p]
      l <- loop_dataframe$j[p]
      feat_k <- k
      feat_l <- l
      cor_k_l <- cor(x[,k], y[,l], method = 'pearson')
      cor_p_val <- cor.test(x[,k], y[,l], method = 'pearson')$p.value
      df <- data.frame(feat_x = feat_k, feat_y = feat_l, cor = cor_k_l, p_val = cor_p_val, stringsAsFactors = FALSE)
      failed <- FALSE
    })
    if (failed) {
      stop(paste('failure with feature pair', feat_k, feat_l))
    }
    return(df)
  }
  print('calculating pairwise correlations')
  loop_df <- data.frame(tidyr::crossing(data.frame(i = colnames(X), stringsAsFactors = FALSE), data.frame(j = colnames(Y), stringsAsFactors = FALSE)),
                        stringsAsFactors = FALSE)

  if (self_only) {
    print('only calculating correlation of features with same name in X and Y')
    self_interaction <- loop_df$i == loop_df$j
    initial_ct <- nrow(loop_df)
    if (sum(self_interaction) == 0) {
      stop('no pairsof features in X and Y features with same name')
    }
    loop_df <- loop_df[self_interaction, ]
    final_ct <- nrow(loop_df)
    print(paste('kept', final_ct, 'of', initial_ct, 'pairs'))
  }
  if (n_cores == 1) {
    df_return <- data.frame(feat_x = character(0), feat_y = character(0), cor = numeric(0), p_val = numeric(0), stringsAsFactors = FALSE)
    pb <- txtProgressBar(max = nrow(loop_df), style = 3)
    ct <- 0
    for (p in 1:nrow(loop_df)) {
      df_return <- rbind(df_return, add_df(X, Y, loop_df, p))
      ct <- ct + 1
      setTxtProgressBar(pb, ct)
    }
    close(pb)
  } else {
    # do parallel runs
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = TRUE)
    bp_result <- BiocParallel::bplapply(1:nrow(loop_df), FUN = function(p) {add_df(X, Y, loop_df, p)})
    df_return <- c()
    for (i in 1:length(bp_result)) {
      df_return <- rbind(df_return, bp_result[[i]])
    }
  }
  cor_fdr <- p.adjust(df_return$p_val, method = 'BH')
  df_return$fdr <- cor_fdr
  df_return <- df_return[df_return$fdr < fdr_filter, ]
  df_return <- df_return[order(df_return$cor, decreasing = TRUE), ]
  return(df_return)
}


#' Correlate Features pulled from SummExpDR object
#' @param x SummExpDR object
#' @param features_x features to pull that will be represented in X in get_cor_vals
#' @param features_y features to pull that will be represented in Y in get_cor_vals
#' @param prune_x regex to remove from values in features_x.
#' use when different data modalities have measurements for same gene and you want to set self_only to TRUE
#' @param prune_y regex to remove from values in features_y.
#' use when different data modalities have measurements for same gene
#' @param assayKey assay to pull if using assay data
#' @param ... other args to get_cor_vals
#' @value data.frane with feature pairs and correlation values
#' @export

setGeneric('correlateFeatures',
           function(x, features_x, features_y, prune_x = NULL, prune_y = NULL, assayKey = NULL, ...) standardGeneric('correlateFeatures'))

setMethod('correlateFeatures',
          signature = 'SummExpDR',
          function(x, features_x, features_y, prune_x = NULL, prune_y = NULL, assayKey = NULL, ...) {
            red_dim_keys <- getReducedDims_keys(x)
            data_fetched <- fetchData(x, c(features_x, features_y), red_dim_keys, assayKey)
            data_x <- as.matrix(data_fetched[,features_x, drop = FALSE])
            data_y <- as.matrix(data_fetched[,features_y, drop = FALSE])
            if (!is.null(prune_x)) {
              colnames(data_x) <- sub(prune_x, '', colnames(data_x))
            }
            if (!is.null(prune_y)) {
              colnames(data_y) <- sub(prune_y, '', colnames(data_y))
            }
            stopifnot(all(is.numeric(data_x)))
            stopifnot(all(is.numeric(data_y)))
            df_return <- get_cor_vals(data_x, data_y, ...)
            return(df_return)
          })
