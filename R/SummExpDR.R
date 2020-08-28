###############################################################################


#' Wrapper around SummarizedExperiment object to contain Reduced Dims
#' as an additional slot
#'
#' Note that subsetting methods are not included as subsetting AFTER performing
#' dimensionality reduction typically means that you won't get the same coordinates
#' if dimensionality reduction is rerun.
#'
#' @slot summ_exp SummarizedExperiment
#' @slot reducedDims list of SummarizedExperiment objects, containing sample names (colData)
#' @export
#' and reduced dim metadata (e.g. loadings, variance explained)

setClass('SummExpDR', representation = list('summ_exp' = 'SummarizedExperiment',
                                            'reducedDims' = 'list'))

# Create Object ---------------------------------------------------------------

create_SummExpDR <- function(summ_exp) {
  # initialize as object with SummarizedExperiment and empty list
  # for reducedDims.
  obj <- new('SummExpDR',
             summ_exp = summ_exp,
             reducedDims = list())
}

# Getters + Setters -----------------------------------------------------------

#' get summarized experiment object from SummExpDR
#' @export

setGeneric('getSummExp', function(x) standardGeneric('getSummExp'))

setMethod('getSummExp',
          'SummExpDR',
          function(x) {
            return(x@summ_exp)
          })

#' get keys for dim reductions
#' @export

setGeneric('getReducedDims_keys', function(x) standardGeneric('getReducedDims_keys'))

setMethod('getReducedDims_keys',
          'SummExpDR',
          function(x) {
            return(names(x@reducedDims))
          })

#' get reduced dims indicated by key
#' @export

setGeneric('getReducedDims', function(x, key) standardGeneric('getReducedDims'))

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
              return(x@reducedDims[[key]])
            }
          })
#' Setter for reducedDims slot
#' @param x SummExpDR
#' @param key character vector of length 1
#' @param value
#' @export

setGeneric('setReducedDims', function(x, key, value) standardGeneric('setReducedDims'))

setMethod('setReducedDims',
          'SummExpDR',
          function(x, key, value) {
            x@reducedDims[[key]] <- value
            return(x)
          })

#' Subset Data
#' @param x SummExpDR object
#' @param rows
#' @param cols
#' @value subsetted SummExpDR object for given rows (features) and cols (samples)
setGeneric('subsetData', function(x, rows, cols) standardGeneric('subsetData'))

setMethod('subsetData',
          'SummExpDR',
          function(x, rows = NULL, cols = NULL) {
            summ_exp <- getSummExp(x)
            dr_keys <- getReducedDims_keys(x)
            all_rows <- rownames(summ_exp)
            all_cols <- colnames(summ_exp)
            if (is.null(rows)) {
              rows <- all_rows
            } else if (is.numeric(rows)) {
              rows <- all_rows[rows]
            } else if (!is.character(rows)) {
              stop('rows must be specified as character or integer')
            }
            if (is.null(cols)) {
              cols <- all_cols
            } else if (is.numeric(cols)) {
              cols <- all_cols[cols]
            } else if (!is.character(cols)) {
              stop('cols must be specified as character or integer')
            }
            summ_exp <- summ_exp[rows, cols]
            new_SummExpDR <- create_SummExpDR(summ_exp)

            for (k in dr_keys) {
              # note that features used to calculate reduced dims, found in
              # rowdata of FactorizedDR objects, are not subsetted
              reduced_dims_k <- getReducedDims(x, k)[,cols]
              new_SummExpDR <- setReducedDims(new_SummExpDR, k, reduced_dims_k)
            }
            return(new_SummExpDR)
          })

# Check Validity --------------------------------------------------------------

.validSummExpDR <- function(x) {
  stopifnot(is(x@summ_exp, 'SummarizedExperiment'))
  stopifnot(is.list(x@reducedDims))
  if (length(x@reducedDims) > 0) {
    for (i in 1:length(x@reducedDims)) {
      stopifnot(is(x@reducedDims[[i]], 'SummarizedExperiment'))
    }
  }
}

setValidity('SummExpDR', .validSummExpDR)
