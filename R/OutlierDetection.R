###############################################################################
### Outlier detection

#' Detect outliers
#'
#' @param object SummExpDR class object, SummarizedExperiment, or mxn numeric
#' data matrix.
#' @param assay.use which assay to use for outlier detection
#' @return if a SummExpDR, a SummExpDR with metadata updated
#' with an outliers column. if a SummarizedExperiment, a SummarizedExperiment
#' with colData updated with an outliers column. if mxn matrix, a logical
#' vector indicating outlier status
#' @export
#'
#' @examples
detect_outliers <- function(object, assay.use = 1) {
  if (is(object, 'matrix') && is.numeric(object)) {
    outlier.result <- rrcovHD::OutlierPCDist(object)
    outlier.status <- logical(nrow(object))
    outlier.status[rrcovHD::getOutliers(outlier.result)] <- TRUE
    return(outlier.status)
  } else if (is(object, 'SummarizedExperiment')) {
    data.use <- t(SummarizedExperiment::assay(object, assay.use))
    outlier.status <- detect_outliers(data.use)
    object$is.outlier <- outlier.status
    return(object)
  } else if (is(object, 'SummExpDR')) {
    object@summ_exp <- detect_outliers(getSummExp(object), assay.use = assay.use)
    return(object)
  } else {
    stop(paste('unrecognized object argument of class', class(object)))
  }
}

#' Filter outliers from data
#'
#' @param object either a SummarizedExperiment object or a SummExpDR class object
#' @description remove outliers from data
#' @return object with outliers removed from data and metadata. if a SummExpDR,
#' all slots except expt.data are destroyed as output is a new SummExpDR with
#' subsetted experimental data and no entries for clustering analyses + differential
#' phenotype analysis
#' @export
#'
#' @examples
remove_outliers <- function(object) {
  if (is(object, 'SummarizedExperiment')) {
    if (!grep('^is.outlier$', colnames(SummarizedExperiment::colData(object)))) {
      stop('outlier detection not run yet for this object')
    }
    if (any(object$is.outlier)) {
      new.expt <- object[, -which(object$is.outlier)]
    }
    return(new.expt)
  } else if (is(object, 'SummExpDR')) {
    object@summ_exp <- remove_outliers(getSummExp(object))
    return(object)
  } else {
    stop(paste('unrecognized object argument of class', class(object)))
  }
}

