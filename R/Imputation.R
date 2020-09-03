###############################################################################
## Data Imputation


#' Wrapper around sklearn.immpute.KNNImputer
#' @param x data, features in rows, samples in columns
#' @param ... other args for sklearn.impute.KNNImputer()
#' @value imputed data matrix
#' @export

KNN_impute <- function(x, ...) {
  imputer <- SKLEARN_IMPUTE$KNNImputer(...)
  x_imputed <- t(imputer$fit_transform(t(x)))
  rownames(x_imputed) <- rownames(x)
  colnames(x_imputed) <- colnames(x)
  return(x_imputed)
}

#' Perform imputation with rescaling of data
#'
#' Data is scaled by multiplying each variable by corresponding entry of scale_factors (scale_factors must)
#'
#' @inheritParams KNN_impute
#' @param scale_factors if NULL, set to 1/stddev for all feats. otherwise, must be numeric
#' of length nrow(x), and names of scale factors must correspond to names of features
#' @value imputed data matrix
#' @export

KNN_scale_impute <- function(x, scale_factors = NULL, ...) {
  feat_means <- matrixStats::rowMeans2(x, na.rm = TRUE)
  names(feat_means) <- rownames(x)
  if (is.null(scale_factors)) {
    scale_factors <- 1/apply(x, MARGIN = 1, FUN = function(x) {sd(x, na.rm = TRUE)})
    names(scale_factors) <- rownames(x)
  }
  stopifnot(length(scale_factors) == nrow(x))
  stopifnot(all(sort(names(scale_factors)) == sort(rownames(x))))
  # center and scale data
  x_hat <- matrix(rep(feat_means[rownames(x)], ncol(x)), byrow = FALSE, nrow = nrow(x))
  scale_mat <- matrix(rep(scale_factors[rownames(x)], ncol(x)), byrow = FALSE, nrow = nrow(x))
  x_scaled <- (x - x_hat)*scale_mat
  x_scaled_imp <- KNN_impute(x_scaled, ...)
  tryCatch({stopifnot(nrow(x_scaled_imp) == nrow(x))},
           error = function(e) {
             stop('number of features in imputed matrix do not match number of features in input matrix, likely have features with no values whatsoever')
           })
  # rescale data to original scale and uncenter data
  x_unscaled_imp <- (x_scaled_imp)/scale_mat + x_hat
  rownames(x_unscaled_imp) <- rownames(x)
  colnames(x_unscaled_imp) <- colnames(x)
  return(x_unscaled_imp)
}
