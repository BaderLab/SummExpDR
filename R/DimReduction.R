###############################################################################
## Dimensionality Reduction Methods


#' Dimensionality Reduction Container
#'
#' inherits from SummarizedExperiment
#' @inherit SummarizedExperiment
#' @slot SummExpDR_assay assay used to create Dimensionality Reduction object.
#' This comes in handy for varimax rotation, where knowledge of which assay
#' used for PCA is necessary to recompute reduced dims.
setClass('DR', contains = 'SummarizedExperiment')

create_DR <- function(assays, ...) {
  summ_exp <- SummarizedExperiment::SummarizedExperiment(assays, ...)
  DR <- new('DR', colData = summ_exp@colData,
            assays = summ_exp@assays,
            NAMES = summ_exp@NAMES,
            elementMetadata = summ_exp@elementMetaData,
            metadata = summ_exp@metadata,
            SummExpDR_assay = NULL)
  return(DR)
}

setMethod('setSummExpDR_assay',
          c('DR'),
          function(x, assay_name) {
            DR@SummExpDR_assay <- assay_name
            return(DR)
          })

#' Get Embeddings
#'
#' @param x DR (Dimensionality Reduction) class object or SummExpDR object
#' @param key key if x is SummExpDR object
#' @param rows = rows to subset (integer or logical indices or character vector)
#' @param cols = cols to subset (integer or logical indices or character vector)
#' @value embeddings matrix

setGeneric('getEmbeddings',  function(x, rows, cols) standardGeneric("getLoadings"))

setMethod('getEmbeddings',
          'DR',
          function(x, rows, cols) {
            return(SummarizedExperiment::assay(x, 1)[rows, cols])
          })

setMethod('getEmbeddings',
          'SummExpDR',
          function(x, key, rows, cols) {
            DR <- getReducedDims(x, key)
            return(getEmbeddings(DR, rows, cols))
          })

#' DR object for Linear Dimensionality Reduction

setClass('LinearDR', contains = 'DR')

create_LinearDR <- function(assays, ...) {
  DR <- create_DR(assays, ...)
  LinearDR <- new('LinearDR', colData = DR@colData,
            assays = DR@assays,
            NAMES = DR@NAMES,
            elementMetadata = DR@elementMetaData,
            metadata = DR@metadata)
  return(LinearDR)
}

#' Get Loadings from Linear Dimensionality Reduction Object
#'
#' @param x LinearDR object or
#' @param red_dims index for reduced dims. Ideally should be character.
#' @param features names of features to use (e.g. gene names), or NULL for all features. If not null should ideally be character
#' @param key key for dimensionality reduction
#' @value Loadings matrix for linear dimensionality reduction

setGeneric('getLoadings',  function(x) standardGeneric("getLoadings"))

setMethod('getLoadings',
          'LinearDR',
          function(x, red_dims = NULL, features = NULL) {
            col_data <- SummarizedExperiment::colData(x)
            # relying on object to not have features with _loading in names of features
            # used for dimensionality reduction
            col_regex <- '_loading$'
            loading_cols <- colnames(col_data)[grep(col_regex, colnames(col_data))]
            if (is.null(red_dims)) {
              red_dims <- rownames(col_data)
            }
            loading_mat <- as.matrix(col_data[red_dims, loading_cols])
            colnames(loading_mat) <- sub(col_regex, '', colnames(loading_mat))
            if (!is.null(features)) {
              loading_mat <- loading_mat[,features]
            }
            return(loading_mat)
          })

setMethod('getLoadings',
          'SummExpDR',
          function(x, key, red_dims = NULL, features = NULL) {
            LinearDR <- getReducedDims(x, key)
            loading_mat <- getLoadings(LinearDR, red_dims, features)
            return(loading_mat)
          })

.LinearDR_validity <- function(x) {
  col_data <- SummarizedExperiment::colData(x)
  all_colnames <- colnames(col_data)
  loading_regex <- '_loading$'
  var_expl_regex <- '^VarExplained$'
  stopifnot(sum(c(grepl(loading_regex), grepl(var_expl_regex))) == ncol(col_data))
}

setValidity('LinearDR', .LinearDR_validity)

# Methods for Running Dimensionality Reduction --------------------------------

#' Run PCA on assay data
#'
#' @param x = SummExpDR, SummarizedExperiment, or matrix.
#' matrix must have the form n features x m samples
#' @param i = assay to pull
#' @param suffix suffix to add, if SummExpDR. By default it's just \'PCA\'
#' @value SummExpDR updated with LinearDR object or LinearDR object otherwise
#'

setGeneric('runPCA',  function(x, i, suffix) standardGeneric("runPCA"))

setMethod('runPCA',
          'matrix',
          function(x, i = NULL, suffix = NULL) {
            # i and suffix are just placeholders
            data_mat <- t(x)
            pca_res <- prcomp(data_mat, center = TRUE, scale. = TRUE)
            pca_coords <- pca_res$x
            loading_mat <- pca_res$rotation
            rownames(loading_mat) <- paste0(rownames(loadings_mat), '_loading')
            pca_var <- (pca_res$sdev)^2
            key <- paste0('pca', suffix)
            SampleID <- rownames(pca_coords)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            row_data <- S4Vectors::DataFrame(t(loading_mat))
            row_data$VarExplained <- pca_var/sum(pca_var)
            pca_summ_exp <- create_LinearDR(assays = list(pca_coords), rowData = row_data, colData = col_data)
            return(pca_summ_exp)
          })

setMethod('runPCA',
          c('SummarizedExperiment'),
          function(x, i, suffix = NULL) {
            pca_summ_exp <- runPCA(SummarizedExperiment::assay(x, i))
            # record assay used
            setSummExpDR_assay(pca_summ_exp, i)
            return(pca_summ_exp)
          })

setMethod('runPCA',
          c('SummExpDR'),
          function(x, i, suffix) {
            pca_summ_exp <- runPCA(getSummExp(x), i)
            key = paste0('PCA', suffix)
            setReducedDims(x, key = key, value = pca_summ_exp)
            return(x)
          })

#' Run varimax rotation on a linear dimensionality reduction
#'
#' @param x SummExpDR object
#' @param key key for LinearDR in ReducedDims to access
#' @param dims_use dims to use (character or integer/logical index, character preferred)
#' @param ... other args to pass to stats::varimax

setGeneric('runVarimax', function(x, key, dims_use = NULL, ...) standardGeneric('runVarimax'))

setMethod('runVarimax',
          'SummExpDR',
          function(x, key, dims_use = NULL, ...) {
            LinearDR <- getReducedDims(x, key)
            loadings_mat <- getLoadings(LinearDR, red_dims = dims_use)
            #loadings_mat is n PCs x p features
            # varimax expects p features x n PCs and returns a loading matrix of that form
            vmax_res <- varimax(t(loadings_mat))
            vmax_loadings <- vmax_res$loadings[,]
            colnames(vmax_loadings) <- paste0('VM', 1:ncol(vmax_loadings))
            # compute new coords
            assay_used <- LinearDR@SummExpDR_assay
            tryCatch(stopifnot(is.character(assay_used)), {
              stop('assay not specified for Linear DR selected')
            })
            summ_exp <- getSummExp(x)
            assay_names <- SummarizedExperiment::assayNames(summ_exp)
            tryCatch(stopifnot(assay_used %in% assay_names), {
              stop('assay not in SummExpDR\'s SummarizedExperiment object\'s list of assay names')
            })
            # coords are in m samples x p features format
            # vmax loadings in p features x n embedding dims format
            data_coords <- t(SummarizedExperiment::assay(summ_exp, assay_used))
            new_coords <- data_coords %*% vmax_loadings
            var_expl <- numeric(ncol(vmax_loadings))
            names(var_expl) <- colnames(vmax_loadings)
            mean_mat <- matrix(rep(colMeans(data_coords), nrow(data_coords)), nrow = nrow(data_coords), byrow = TRUE)
            total_error <- sum(new_)
            # for (v_name in names(var_expl)) {
            #   reconstr_data <- new_coords[,v_name, drop = FALSE] %*% t(vmax_loadings[,v_name, drop = FALSE])
            #   diff <- reconstr_data - data_coords
            #   sq_error <- sum(diag(diff %*% t(diff))^2)
            # }
            # reconstr_data <- new_coords[,2:3, drop = FALSE] %*% t(vmax_loadings[,2:3])
            # diff <- reconstr_data - data_coords
            # sq_error <- sum(diag(diff %*% t(diff))^2)

          })
