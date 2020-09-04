###############################################################################
## Multi-assay integration

#' object for integrated analysis of multiple datatypes
#' @include SummExpDR.R
#' @inherit SummExpDR
#' @export

setClass('multiExp',
         contains = 'SummExpDR')

# create multiExp -----------------------------------------------------------

#' create multiExp object
#'
#' Creates a new SummarizedExperiment object for summ_exp slot, containing combined data from multiple experiments
#' E.g. RNA, DNA methylation, CNV, etc. Data is standard normalized prior to performing imputation. Choice of
#' word 'Experiment' to describe multiple data modalities is to avoid confusion with concept of 'assays' in
#' SummarizedExperiment, which assumes identical set of biological feature names with different measurements
#'
#' @param summ_exp_list list of SummarizedExperiment objects
#' @param assays_use character vector or numeric vector with names corresponding to names in summ_exp_list
#' @param ... args to pass to KNN_impute()
#' @value multiExp object
#' @export

createMultiExp <- function(summ_exp_list, assays_use = NULL, ...) {



  expt_names <- names(summ_exp_list)

  tryCatch({
    stopifnot(is.character(expt_names))
  }, error = function(e) {
    stop('names for experiments must be specified in summ_exp_list')
  })

  tryCatch({
    stopifnot(all(unlist(lapply(summ_exp_list, FUN = function(x) {is(x, 'SummarizedExperiment')}))))
  }, error = function(e) {
   stop('all summ_exp_list elements mut be of class SummarizedExperiment')
  })

  if (is.null(assays_use)) {
    warning('assays_use not specified. using 1st assay of each provided SummarizedExperiment')
    assays_use <- rep(1L, length(summ_exp_list))
    names(assays_use) <- expt_names
  } else {
    tryCatch({
      stopifnot(is.character(names(assays_use)))
      stopifnot(sort(names(assays_use)) == sort(expt_names))
    }, error = function(e) {
      stop('names for assays_use must be specified and must match summ_exp_list names')
    })
  }

 feature_key <- S4Vectors::DataFrame(expt = character(0), orig_id = character(0), converted = character(0))

  for (i in 1:length(summ_exp_list)) {
    summ_exp_i <- check_rownames_colnames(summ_exp_list[[i]])
    # rename features
    expt_name_i <- expt_names[i]
    orig_names <- rownames(summ_exp_i)
    new_names <- paste(expt_name_i, orig_names, sep = '_')
    rownames(summ_exp_i) <- new_names
    # save name mapping between original feature names and new feature names.
    # this is added to row data
    feature_key_i <- S4Vectors::DataFrame(expt = rep(expt_name_i, nrow(summ_exp_i)),
                                          orig_id = orig_names,
                                          converted = new_names)
    rownames(feature_key_i) <- new_names
    assay_i <- assays_use[[i]]
    assay_mat_i <- SummarizedExperiment::assay(summ_exp_i, assay_i)
    # get row and coldata
    col_data_i <- SummarizedExperiment::colData(summ_exp_i)
    row_data_i <- SummarizedExperiment::rowData(summ_exp_i)
    row_data_i <- replace_col(row_data_i, 'feat_id', rownames(row_data_i), suffix = paste0('_', expt_name_i))
    col_data_i <- replace_col(col_data_i, 'sample_id', rownames(col_data_i), suffix = paste0('_', expt_name_i))

    # row_data_i <- cbind(feature_key_i, row_data_i)
    feature_key <- rbind(feature_key, feature_key_i)
    if (i == 1) {
      new_mat <- assay_mat_i
      row_data <- row_data_i
      col_data <- col_data_i
    } else {
      # there cannot be intersecting rownames since they are prefixed with
      # experiment name, which is unique. thus no duplicate experiment/gene combinations
      # we do have to handle missing samples from current new_mat and current assay_mat_i
      # (i.e. non intersecting samples), and we do this by inserting NA values
      missing_cols_x <- setdiff(colnames(new_mat), colnames(assay_mat_i))
      missing_cols_y <- setdiff(colnames(assay_mat_i), colnames(new_mat))
      add_na_vals <- function(mat, col_names) {
        mat_append <- matrix(data = NA, nrow = nrow(mat), ncol = length(col_names))
        rownames(mat_append) <- rownames(mat)
        colnames(mat_append) <- col_names
        mat <- cbind(mat, mat_append)
      }
      if (length(missing_cols_x) > 0) {
        assay_mat_i <- add_na_vals(assay_mat_i, missing_cols_x)
      }
      if (length(missing_cols_y) > 0) {
        new_mat <- add_na_vals(new_mat, missing_cols_y)
      }
      # join data matrix, rowdata, coldata
      new_mat <- rbind(new_mat[ , sort(colnames(new_mat))], assay_mat_i[ , sort(colnames(new_mat))])
      row_data <- S4Vectors::DataFrame(dplyr::full_join(as.data.frame(row_data),
                                                        as.data.frame(row_data_i),
                                                        by = 'feat_id',
                                                        suffix = c('', paste0('_', expt_name_i))
                                                        )
                                       )
      rownames(row_data) <- row_data[['feat_id']]
      col_data <- S4Vectors::DataFrame(dplyr::full_join(as.data.frame(col_data),
                                                        as.data.frame(col_data_i),
                                                        by = 'sample_id',
                                                        suffix = c('', paste0('_', expt_name_i))
                                                        )
                                       )
      rownames(col_data) <- col_data[['sample_id']]
    }
    # num_features[expt_name_i] <- nrow(assay_mat_i)
  }
  rownames(feature_key) <- feature_key$converted
  for (col_name in c('expt', 'orig_id')) {
    row_data <- replace_col(row_data, col_name, as.vector(feature_key[rownames(row_data), col_name]), suffix = paste0('_orig'))
  }
  new_mat <- new_mat[ , sort(colnames(new_mat))]
  col_data <- col_data[sort(rownames(col_data)), ]
  new_summ_exp <- SummarizedExperiment::SummarizedExperiment(assays = list(stacked = new_mat), rowData = row_data, colData = col_data)
  newSummExpDR <- create_SummExpDR(new_summ_exp)
  newMultiExp <- new('multiExp',
                     summ_exp = newSummExpDR@summ_exp,
                     reducedDims = newSummExpDR@reducedDims)
  return(newMultiExp)
}

setGeneric('getScaleFactors',  function(x, assay) standardGeneric("getScaleFactors"))

setMethod('getScaleFactors',
          signature = 'multiExp',
          function(x, assay) {
            summ_expt <- getSummExp(x)
            data_mat <- SummarizedExperiment::assay(summ_expt, assay)
            row_data <- SummarizedExperiment::rowData(summ_expt)
            unique_expts <- unique(row_data$expt)
            scale_factors <- c()
            for (expt_name in unique_expts) {
              feat_id <- row_data[row_data$expt == expt_name, 'feat_id']
              scale_factors_expt <- numeric(length(feat_id))
              names(scale_factors_expt) <- feat_id
              scale_factors_expt[1:length(scale_factors_expt)] <- apply(data_mat[feat_id,],
                                                                        MARGIN = 1,
                                                                        FUN = function(x) {sd(x, na.rm = TRUE)})
              # scale factor = 1/(sd(feat)*sqrt(# feats in experiment))
              scale_factors_expt <- 1/(scale_factors_expt*sqrt(length(feat_id)))
              scale_factors <- c(scale_factors, scale_factors_expt)
            }
            return(scale_factors)
          })

#' Impute Missing Data in multiExp object
#' @param x multiExp object
#' @param assay assay for which data is to be imputed. by default it's stacked, the assay name for all of the input data stacked together
#' @param ... other args for sklearn.impute.KNNImputer
#' @export

setGeneric('imputeExpData',  function(x, assay = 'stacked', ...) standardGeneric("imputeExpData"))

setMethod('imputeExpData',
          signature = 'multiExp',
          function(x, assay = 'stacked', ...) {
            summ_expt <- getSummExp(x)
            data_mat <- SummarizedExperiment::assay(summ_expt, assay)
            na_inds <- which(is.na(data_mat))
            is_imputed <- matrix(data = 0L, nrow = nrow(data_mat), ncol = ncol(data_mat))
            rownames(is_imputed) <- rownames(data_mat)
            colnames(is_imputed) <- colnames(data_mat)
            if (length(na_inds) > 0) {
              print(paste('running imputation for', length(na_inds), 'of', length(data_mat), 'values'))
              is_imputed[na_inds] <- 1L
              scale_factors <- getScaleFactors(x, assay)
              imputed_mat <- KNN_scale_impute(data_mat, scale_factors, ...)
              SummarizedExperiment::assay(summ_expt, 'imputed_mat') <- imputed_mat
            } else {
              print('nothing to impute')
            }
            SummarizedExperiment::assay(summ_expt, 'is_imputed') <- is_imputed
            x@summ_exp <- summ_expt
            return(x)
          })


#' Scale Experiment Data
#'
#' @param x multiExp object
#' @param assay character, assay to pull out for scaling.
#' @value multiExp object with 'scaled' data added
#' @export

setGeneric('scaleExpData', function(x, assay) standardGeneric('scaleExpData'))

setMethod('scaleExpData',
          signature = 'multiExp',
          function(x, assay) {
            # standard deviation is same whether it's computed before or after centering data
            # recall that scale factor is 1/(sd * sqrt(feats))
            scale_factors <- getScaleFactors(x, assay)
            summ_expt <- getSummExp(x)
            data_mat <- SummarizedExperiment::assay(summ_expt, assay)
            data_mat <- t(scale(t(data_mat), scale = FALSE, center = TRUE))
            scale_mat <- matrix(rep(scale_factors[rownames(data_mat)], ncol(data_mat)),
                                nrow = nrow(data_mat),
                                byrow = FALSE)
            scaled_data <- data_mat*scale_mat
            SummarizedExperiment::assay(x@summ_exp, 'scaled') <- scaled_data
            return(x)
          })
