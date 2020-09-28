###############################################################################


#' Wrapper around SummarizedExperiment object to contain Reduced Dims
#' as an additional slot
#'
#' Note that subsetting methods are not included as subsetting AFTER performing
#' dimensionality reduction typically means that you won't get the same coordinates
#' if dimensionality reduction is rerun.
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @slot summ_exp SummarizedExperiment
#' @slot reducedDims list of SummarizedExperiment objects, containing sample names (colData)
#' and reduced dim metadata (e.g. loadings, variance explained)
#' @export


setClass('SummExpDR', representation = list('summ_exp' = 'SummarizedExperiment',
                                            'reducedDims' = 'list'))

# Create Object ---------------------------------------------------------------

check_rownames_colnames <- function(x) {
  stopifnot(is(x, 'SummarizedExperiment'))
  orig_rownames <- rownames(x)
  orig_colnames <- colnames(x)
  rownames(x) <- make.names(orig_rownames)
  colnames(x) <- make.names(orig_colnames)
  if (!all(rownames(x) == orig_rownames)) {
    warning('some or all rownames in input were changed as they had spaces or disallowed symbols')
  }
  row_replace <- orig_rownames
  SummarizedExperiment::rowData(x) <- replace_col(SummarizedExperiment::rowData(x),
                                                  col_name = 'raw_rownames', value = row_replace, suffix = 'orig')
  if (!all(colnames(x) == orig_colnames)) {
    warning('some or all colnames in input were changed as they had spaces or disallowed symbols')
  }
  col_replace <- orig_colnames
  SummarizedExperiment::colData(x) <- replace_col(SummarizedExperiment::colData(x),
                                                  col_name = 'raw_colnames', value = col_replace, suffix = 'orig')
  return(x)
}

#' Create SummExpDR object
#' @param summ_exp SummarizedExperiment object
#' @export

create_SummExpDR <- function(summ_exp) {
  # initialize as object with SummarizedExperiment and empty list
  # for reducedDims.
  summ_exp <- check_rownames_colnames(summ_exp)
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

#' Row Data
#' @param x SummExpDR object
#' @value rowData from Summ Exp DR object
#' @export
setGeneric('rowData', function(x) standardGeneric('rowData'))

setMethod('rowData',
          signature = 'SummExpDR',
          function(x) {
            return(SummarizedExperiment::rowData(getSummExp(x)))
          })

#' add rowData
#' @param x
#' @param value DataFrame or coercible thereof, or vector. if vector, col_name must be specified
#' @param col_name
#' @value SummExpDR object with rowData updated
#' @export

setGeneric('addRowData', function(x, value, col_name = NULL) standardGeneric('addRowData'))

setMethod('addRowData',
          signature = 'SummExpDR',
          function(x, value, col_name = NULL) {
            if (is.vector(value)) {
              if (is.null(col_name)) {
                stop('col_name must not be null if value is a vector')
              }
              if (is.null(names(value))) {
                stop('vector input for value must have bames')
              }
              value <- S4Vectors::DataFrame(x = value, row.names = names(value))
              colnames(value) <- col_name
            }
            row_data <- SummarizedExperiment::rowData(x@summ_exp)
            missing_names_rowdata <- setdiff(rownames(row_data), rownames(value))
            missing_names_value <- setdiff(rownames(value), rownames(row_data))
            if (length(missing_names_rowdata) > 0) {
              warning(paste(length(missing_names_rowdata) , 'rownames in rowdata not represented in value argument'))
            }
            if (length(missing_names_value) > 0) {
              warning(paste(length(missing_names_value) , 'rownames/names in value not represented in rowdata'))
            }
            for (v in colnames(value)) {
              row_data <- replace_col(row_data, v, value[rownames(row_data),v], suffix = 'orig', remove_existing = TRUE)
            }
            SummarizedExperiment::rowData(x@summ_exp) <- col_data
            return(x)
          })

#' Col Data
#' @param x SummExpDR object
#' @value rowData from Summ Exp DR object
#' @export
setGeneric('colData', function(x) standardGeneric('colData'))

setMethod('colData',
          signature = 'SummExpDR',
          function(x) {
            return(SummarizedExperiment::colData(getSummExp(x)))
          })


#' add colData
#' @param x
#' @param value DataFrame or coercible thereof, or vector. if vector, col_name must be specified
#' @param col_name
#' @value SummExpDR object with colData updated
#' @export

setGeneric('addColData', function(x, value, col_name = NULL) standardGeneric('addColData'))

setMethod('addColData',
          signature = 'SummExpDR',
          function(x, value, col_name = NULL) {
            if (is.vector(value)) {
              if (is.null(col_name)) {
                stop('col_name must not be null if value is a vector')
              }
              if (is.null(names(value))) {
                stop('vector input for value must have bames')
              }
              value <- S4Vectors::DataFrame(x = value, row.names = names(value))
              colnames(value) <- col_name
            }
            col_data <- SummarizedExperiment::colData(x@summ_exp)
            missing_names_coldata <- setdiff(rownames(col_data), rownames(value))
            missing_names_value <- setdiff(rownames(value), rownames(col_data))
            if (length(missing_names_coldata) > 0) {
              warning(paste(length(missing_names_coldata) , 'rownames in coldata not represented in value argument'))
            }
            if (length(missing_names_value) > 0) {
              warning(paste(length(missing_names_value) , 'rownames/names in value not represented in coldata'))
            }
            for (v in colnames(value)) {
              col_data <- replace_col(col_data, v, value[rownames(col_data),v], suffix = 'orig', remove_existing = TRUE)
            }
            SummarizedExperiment::colData(x@summ_exp) <- col_data
            return(x)
          })

#' Subset Data
#' @param x SummExpDR object
#' @param rows
#' @param cols
#' @value subsetted SummExpDR object for given rows (features) and cols (samples)
#' @export
setGeneric('subsetData', function(x, rows = NULL, cols = NULL) standardGeneric('subsetData'))

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
            x@summ_exp <- summ_exp

            for (k in dr_keys) {
              # note that features used to calculate reduced dims, found in
              # rowdata of FactorizedDR objects, are not subsetted
              reduced_dims_k <- getReducedDims(x, k)[,cols]
              x <- setReducedDims(x, k, reduced_dims_k)
            }
            return(x)
          })

#' Get Assays
#' @param x SummExpDR object
#' @value list of assays
#' @export
setGeneric('assays', function(x) standardGeneric('assays'))

setMethod('assays',
          signature = 'SummExpDR',
          function(x) {
            return(SummarizedExperiment::assays(getSummExp(x)))
          })

#' Show Assay Names
#' @param x SummExpDR object
#' @value vector of assay names
#' @export
setGeneric('assayNames', function(x) standardGeneric('assayNames'))

setMethod('assayNames',
          signature = 'SummExpDR',
          function(x) {
            return(SummarizedExperiment::assayNames(getSummExp(x)))
          })

#' Get Assay Data
#' @param x SummExpDR object
#' @param i assay data to pull
#' @value assay data
#' @export
setGeneric('assay', function(x, i) standardGeneric('assay'))

setMethod('assay',
          signature = 'SummExpDR',
          function(x, i) {
            return(SummarizedExperiment::assay(getSummExp(x), i))
          })

#' helper function to replace column
replace_col <- function(DF, col_name, value, suffix, remove_existing = FALSE) {
  if (col_name %in% colnames(DF)) {
    if (remove_existing) {
      warning(paste(col_name, 'in column names of DF replaced with new value'))
    } else {
      new_colname <- paste(col_name, suffix, sep = '_')
      colnames(DF)[grep(paste0('^', col_name, '$'), col_name)] <- new_colname
      warning(paste(col_name, 'in column names of DF renamed to', new_colname))
    }
  }
  DF[[col_name]] <- value
  return(DF)
}

#' Interface to fetch data from object
#' @param x SummExpDR object
#' @param varsFetch variables to pull out
#' @param redDimKeys keys for reduced dims to pull out
#' @param assayKey assay to pull from summExp object
#' @param mode whether to pull data in a sample x feature format or feature x feature_info format
#' for sample x feature (sample_wises), first pull colData of summ exp object, then join with embeddings from redDimKeys
#' and with assay data from assay keys. for feature x feat info, pull out rowData, then any loadings for reduced dim keys
#' @value data.frame with fetched data
#' @export

setGeneric('fetchData', function(x, varsFetch, redDimKeys = NULL, assayKey = NULL, mode = 'sample_wise') standardGeneric('fetchData'))

setMethod('fetchData',
          signature = 'SummExpDR',
          function(x, varsFetch, redDimKeys = NULL, assayKey = NULL, mode = 'sample_wise') {
            summ_exp <- getSummExp(x)
            red_dim_list <- list()
            if (!is.null(redDimKeys)) {
              for (k in redDimKeys) {
                red_dim_list[[k]] <- getReducedDims(x, k)
              }
            }

            if (mode == 'sample_wise') {
              data_fetched <- as.data.frame(SummarizedExperiment::colData(summ_exp))
              data_fetched <- replace_col(data_fetched, col_name = 'row_names',
                                          value = rownames(data_fetched), suffix = 'colData')

              if (length(red_dim_list) > 0) {
                for (k in redDimKeys) {
                  embeddings_k <- as.data.frame(t(getEmbeddings(red_dim_list[[k]])))
                  embeddings_k <- replace_col(embeddings_k, col_name = 'row_names',
                                              value = rownames(embeddings_k), suffix = k)
                  # tryCatch({stopifnot(any(varsFetch %in% colnames(embeddings_k)))},
                  #          error = function(e) {
                  #            stop('no variables specified found in selected reduced dims\' embeddings')
                  #          })
                  data_fetched <- dplyr::left_join(data_fetched, embeddings_k, by = 'row_names', suffix = c('', k))
                }
              }

              if (!is.null(assayKey)) {
                feat_names <- rownames(summ_exp)
                feats_extract <- intersect(feat_names, varsFetch)
                # tryCatch({stopifnot(length(feats_extract) > 0)},
                #          error = function(e) {
                #            stop('no variables specified found in assay feature names')
                #          })
                assay_data <- as.data.frame(t(SummarizedExperiment::assay(summ_exp, assayKey)))[,feats_extract, drop = FALSE]
                assay_data <- replace_col(assay_data, col_name = 'row_names',
                                            value = rownames(assay_data), suffix = assayKey)
                data_fetched <- dplyr::left_join(data_fetched, assay_data, by = 'row_names', )
              }

            } else if (mode == 'feature_wise') {
              data_fetched <- as.data.frame(SummarizedExperiment::rowData(summ_exp))
              data_fetched <- replace_col(data_fetched, col_name = 'row_names',
                                          value = rownames(data_fetched), suffix = 'rowData')
              if (length(red_dim_list) > 0) {
                for (k in redDimKeys) {
                  loadings_k <- as.data.frame(t(getLoadings(red_dim_list[[k]])))
                  loadings_k <- replace_col(loadings_k, col_name = 'row_names',
                                              value = rownames(loadings_k), suffix = k)
                  # tryCatch({stopifnot(any(varsFetch %in% colnames(loadings_k)))},
                  #          error = function(e) {
                  #            stop('no variables specified found in selected reduced dims\' embeddings')
                  #          })
                  data_fetched <- dplyr::left_join(data_fetched, loadings_k, by = 'row_names', suffix = c('', k))
                }
              }
            } else {
              stop('mode must be one of sample_wise, feature_wise')
            }
            tryCatch({stopifnot(all(varsFetch %in% colnames(data_fetched)))},
                     error = function(e) {
                       stop(paste('following columns missing from fetched data:',
                                  paste(setdiff(varsFetch, colnames(data_fetched)), collapse = ',')))
                     })
            rownames(data_fetched) <- data_fetched[['row_names']]
            data_fetched <- data_fetched[ ,varsFetch, drop = FALSE]
            return(data_fetched)
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
