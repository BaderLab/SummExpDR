###############################################################################


#' Wrapper around SummarizedExperiment object to contain Reduced Dims
#' as an additional slot
#'
#' Note that subsetting methods are not included as subsetting AFTER performing
#' dimensionality reduction typically means that you won't get the same coordinates
#' if dimensionality reduction is rerun.
#'
#' @slot summ_exp SummarizedExperiment
#' @slot ReducedDims list of SummarizedExperiment objects, containing sample names (colData)
#' and reduced dim metadata (e.g. loadings, variance explained)

setClass('SummExpDR', representation = c('summ_exp' = 'SummarizedExperiment',
                                         'ReducedDims' = 'list'))

# Create Object ---------------------------------------------------------------

create_SummExpDR <- function(summ_exp) {
  # initialize as object with SummarizedExperiment and empty list
  # for ReducedDims.
  obj <- new('SummExpDR',
             summ_exp = summ_exp,
             ReducedDims = list())
}

# Getters + Setters -----------------------------------------------------------

setMethod('getSummExp',
          'SummExpDR',
          function(x) {
            return(x@summ_exp)
          })

setMethod('getReducedDims_keys',
          'SummExpDR',
          function(x) {
            return(names(x@reducedDims))
          })

setMethod('getReducedDims',
          'SummExpDR',
          function(x, key = NULL) {
            if (!(is.character(key) && length(key) == 1)) {
              stop('key must be character of length 1')
            } else {
              valid_keys <- getReducedDims_keys(x)
              if (!key %in% valid_keys) {
                stop(paste('provided key not in keys', paste(valid_keys, collapse = ',')))
              }
              return(x@ReducedDims[[key]])
            }
          })
#' Setter for ReducedDims slot
#' @param x SummExpDR
#' @param key character vector of length 1
#' @param value
setMethod('setReducedDims',
          'SummExpDR',
          function(x, key, value) {
            x@ReducedDims[[key]] <- value
          })

# Check Validity --------------------------------------------------------------

.validSummExpDR <- function(x) {
  stopifnot(is(x@summ_exp, 'SummarizedExperiment'))
  stopifnot(is.list(x@ReducedDims))
  if (length(x@ReducedDims) > 0) {
    for (i in 1:length(x@ReducedDims)) {
      stopifnot(is(x@ReducedDims[[i]], 'SummarizedExperiment'))
    }
  }
}

setValidity('SummExpDR', .validSummExpDR)
