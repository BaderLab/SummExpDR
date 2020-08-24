###############################################################################
## Dimensionality Reduction Methods


#' Dimensionality Reduction Container
#'
#' inherits from SummarizedExperiment. assay data = dimensionality reduction
#' coordinates, rows correspond to PCs, cols correspond to samples
#' @inherit SummarizedExperiment::SummarizedExperiment
#' @slot sourceAssay assay used to create Dimensionality Reduction object.
#' This comes in handy for varimax rotation, where knowledge of which assay
#' used for PCA is necessary to recompute reduced dims.
#' @export
setClass('DR', contains = 'SummarizedExperiment', representation(sourceAssay = 'character_OR_NULL'))

create_DR <- function(assays, ...) {
  summ_exp <- SummarizedExperiment::SummarizedExperiment(assays, ...)
  DR <- new('DR', colData = summ_exp@colData,
            assays = summ_exp@assays,
            NAMES = summ_exp@NAMES,
            elementMetadata = summ_exp@elementMetadata,
            metadata = summ_exp@metadata,
            sourceAssay = NULL)
  return(DR)
}

#' Set assay assay used to generate dim reduction in SummExpDR object
#' @export

setGeneric('setSourceAssay', function(x, assay_name) standardGeneric('setSourceAssay'))

setMethod('setSourceAssay',
          c('DR'),
          function(x, assay_name) {
            DR <- x
            DR@sourceAssay <- assay_name
            return(DR)
          })

#' Get Embeddings
#'
#' @param x DR (Dimensionality Reduction) class object or SummExpDR object
#' @param key key if x is SummExpDR object
#' @param rows = rows to subset (integer or logical indices or character vector)
#' @param cols = cols to subset (integer or logical indices or character vector)
#' @value embeddings matrix
#' @export

setGeneric('getEmbeddings',  function(x, key = NULL, rows = NULL, cols = NULL) standardGeneric("getEmbeddings"))

setMethod('getEmbeddings',
          'DR',
          function(x, key = NULL, rows = NULL, cols = NULL) {
            if (is.null(rows)) {
              rows <- 1:nrow(x)
            }
            if (is.null(cols)) {
              cols <- 1:ncol(x)
            }
            return(SummarizedExperiment::assay(x, 1)[rows, cols])
          })

setMethod('getEmbeddings',
          'SummExpDR',
          function(x, key, rows = NULL, cols = NULL) {
            DR <- getReducedDims(x, key)
            return(getEmbeddings(DR, rows = rows, cols = cols))
          })

#' DR object for Linear Dimensionality Reduction
#' @export
setClass('LinearDR', contains = 'DR')

#' Create LinearDR object
#'@export

create_LinearDR <- function(assays, ...) {
  DR <- create_DR(assays, ...)
  LinearDR <- new('LinearDR', colData = DR@colData,
            assays = DR@assays,
            NAMES = DR@NAMES,
            elementMetadata = DR@elementMetadata,
            metadata = DR@metadata,
            sourceAssay = NULL)
  return(LinearDR)
}

#' Get Loadings from Linear Dimensionality Reduction Object
#'
#' @param x LinearDR object or
#' @param red_dims index for reduced dims. Ideally should be character.
#' @param features names of features to use (e.g. gene names), or NULL for all features. If not null should ideally be character
#' @param key key for dimensionality reduction
#' @value Loadings matrix for linear dimensionality reduction
#' @export

setGeneric('getLoadings',  function(x, key = NULL, red_dims = NULL, features = NULL) standardGeneric("getLoadings"))

setMethod('getLoadings',
          'LinearDR',
          function(x, key = NULL, red_dims = NULL, features = NULL) {
            row_data <- SummarizedExperiment::rowData(x)
            # relying on object to not have features with _loading in names of features
            # used for dimensionality reduction
            col_regex <- '_loading$'
            loading_cols <- colnames(row_data)[grep(col_regex, colnames(row_data))]
            if (is.null(red_dims)) {
              red_dims <- 1:nrow(row_data)
            }
            loadings_mat <- as.matrix(row_data[red_dims, loading_cols])
            colnames(loadings_mat) <- sub(col_regex, '', colnames(loadings_mat))
            if (!is.null(features)) {
              loadings_mat <- loadings_mat[,features]
            }
            return(loadings_mat)
          })

setMethod('getLoadings',
          'SummExpDR',
          function(x, key, red_dims = NULL, features = NULL) {
            LinearDR <- getReducedDims(x, key)
            key <- NULL
            loadings_mat <- getLoadings(LinearDR, red_dims = red_dims, features = features)
            return(loadings_mat)
          })

.LinearDR_validity <- function(x) {
  row_data <- SummarizedExperiment::rowData(x)
  all_colnames <- colnames(row_data)
  loading_regex <- '_loading$'
  stopifnot(sum(grepl(loading_regex, all_colnames)) == ncol(row_data))
}

setValidity('LinearDR', .LinearDR_validity)

#' class for linear dimensionality reduction that
#' can have data reconstructed from factors
#' @export
setClass('FactorizedDR', contains = 'LinearDR')

create_FactorizedDR <- function(assays, ...) {
  LinearDR <- create_LinearDR(assays, ...)
  FactorizedDR <- new('FactorizedDR', colData = LinearDR@colData,
                  assays = LinearDR@assays,
                  NAMES = LinearDR@NAMES,
                  elementMetadata = LinearDR@elementMetadata,
                  metadata = LinearDR@metadata,
                  sourceAssay = NULL)
  return(FactorizedDR)
}

#' #' method for calculating reconstruction accuracy
#' #'
#'
#' setGeneric('qualityMetric',  function(x, key, dims_use) standardGeneric("qualityMetric"))
#'
#' setMethod('qualityMetric',
#'           signature = 'SummExpDR',
#'           function(x, key, dims_use) {
#'             FactorizedDR <- getReducedDims(x, key)
#'             assay_used <- FactorizedDR@sourceAssay
#'             input_data_SE <- getSummExp()
#'           })

# Methods for Running Dimensionality Reduction --------------------------------

#' Run PCA on assay data
#'
#' @param x = SummExpDR, SummarizedExperiment, or matrix.
#' matrix must have the form n features x m samples
#' @param i = assay to pull
#' @param suffix suffix to add, if SummExpDR. By default it's just \'PCA\'
#' @param std_norm whether to standard normalize data
#' @value SummExpDR updated with LinearDR object or LinearDR object otherwise
#' @export

setGeneric('runPCA',  function(x, i = NULL, suffix = NULL, std_norm = TRUE) standardGeneric("runPCA"))

setMethod('runPCA',
          'matrix',
          function(x, i = NULL, suffix = NULL, std_norm = TRUE) {
            # i and suffix are just placeholders
            data_mat <- t(x)
            if (std_norm) {
              data_mat <- scale(data_mat, center = TRUE, scale = TRUE)
            }
            # tryCatch({
            #   stopifnot(all(Matrix::colMeans(x) == 0))
            #   stopifnot(all(matrixStats::colVars(x) == 1))
            # }, error = function(e) {
            #   stop('features not standard normalized')
            # })
            pca_res <- prcomp(data_mat, center = FALSE, scale. = FALSE)
            pca_coords <- pca_res$x
            loadings_mat <- pca_res$rotation
            rownames(loadings_mat) <- paste0(rownames(loadings_mat), '_loading')
            # pca_var <- (pca_res$sdev)^2
            key <- paste0('pca', suffix)
            SampleID <- rownames(pca_coords)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            row_data <- S4Vectors::DataFrame(t(loadings_mat))
            # row_data$VarExplained <- pca_var/sum(pca_var)
            # recall that SummarizedExperiment objects have feats in rows and samples in cols
            pca_summ_exp <- create_FactorizedDR(assays = list(pca = t(pca_coords)), rowData = row_data, colData = col_data)
            return(pca_summ_exp)
          })

setMethod('runPCA',
          c('SummarizedExperiment'),
          function(x, i, suffix = NULL, std_norm = TRUE) {
            if (is.numeric(i)) {
              assay_name <- SummarizedExperiment::assays(x)[i]
            } else {
              assay_name <- i
            }
            if (std_norm) {
              new_assay_name <- paste0(assay_name, '_std_norm')
              SummarizedExperiment::assay(x, new_assay_name) <- t(scale(t(SummarizedExperiment::assay(x, assay_name))))
              assay_use <- new_assay_name
            } else {
              assay_use <- assay_name
            }
            pca_summ_exp <- runPCA(SummarizedExperiment::assay(x, assay_use), std_norm = FALSE)
            # record assay used
            setSourceAssay(pca_summ_exp, assay_use)
            return(pca_summ_exp)
          })

setMethod('runPCA',
          c('SummExpDR'),
          function(x, i, suffix, std_norm = TRUE) {
            pca_summ_exp <- runPCA(getSummExp(x), i, suffix, std_norm)
            key = paste0('PCA', suffix)
            x <- setReducedDims(x, key = key, value = pca_summ_exp)
            return(x)
          })

#' Run varimax rotation on a linear dimensionality reduction
#'
#' @param x SummExpDR object
#' @param key key for LinearDR in ReducedDims to access
#' @param dims_use dims to use (character or integer/logical index, character preferred)
#' @param ... other args to pass to stats::varimax
#' @export

setGeneric('runVarimax', function(x, key, dims_use = NULL, suffix = '', ...) standardGeneric('runVarimax'))

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
            assay_used <- LinearDR@sourceAssay
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
            # we transpose these in constructing new LinearDR (SummarizedExperiment) object
            data_coords <- t(SummarizedExperiment::assay(summ_exp, assay_used))
            new_coords <- data_coords %*% vmax_loadings
            rownames(new_coords) <- rownames(data_coords)
            row_data = S4Vectors::DataFrame(as.data.frame(t(loadings)))
            SampleID <- rownames(new_coords)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            VmaxDR <- create_LinearDR(assays = list(varimax = t(new_coords)),
                                      rowData = row_data,
                                      colData = col_data)
            vmax_key <- paste0('varimax', suffix)
            x <- setReducedDims(x, key = vmax_key, value = VmaxDR)
            return(x)
          })
