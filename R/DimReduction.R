###############################################################################
## Dimensionality Reduction Methods


#' Dimensionality Reduction Container
#'
#' inherits from SummarizedExperiment. assay data = dimensionality reduction
#' coordinates, rows correspond to PCs, cols correspond to samples
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @include SummExpDR.R
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

#' method for calculating reconstruction accuracy
#'
#' Calculates reconstruction accuracy (r-squared) of matrix factorization w.r.t. input data for
#' fitting matrix factorization
#'
#' @param x SummExpDR with FactorizedDR dimensionality reduction object computed
#' @param key string indicating which item from reducedDims to compute reconstruction accuracy for
#' @param dims_use integer vector or character vector indicating dimensions to use, or NULL for all dimensions
#' @param feats_use integer vector or character vector indicating features to compute variance explained for. if NULL,
#' all features. Essentially you subset the data matrix and the reconstructed data matrix for the selected features,
#' then compute r-squared to compare those subset matrices.
#' @value returns list with variance explained per dimension as well as total variance explained by dimensions selected
#' @export

setGeneric('varianceExplained',  function(x, key, dims_use = NULL, feats_use = NULL) standardGeneric("varianceExplained"))

setMethod('varianceExplained',
          signature = 'SummExpDR',
          function(x, key, dims_use = NULL, feats_use = NULL) {
            FactorizedDR <- getReducedDims(x, key)
            assay_used <- FactorizedDR@sourceAssay
            # get input data used for this dim reduction
            input_data_SE <- getSummExp(x)
            input_data <- t(SummarizedExperiment::assay(input_data_SE, assay_used))

            # loadings matrix in in p dimensions x n features format
            loadings_mat <- getLoadings(FactorizedDR)
            tryCatch(stopifnot(all(colnames(loadings_mat) == colnames(input_data))),
                     error = function(e) {
                       stop('colnames of loadings matrix do not match feature names of input data')
                     })


            # get embeddings in m samples x p dimensions format
            embeddings_mat <- t(getEmbeddings(FactorizedDR))
            # setup variable to contain variance by dimension
            dim_names <- rownames(loadings_mat)
            if (is.null(dims_use)) {
              dims_use <- dim_names
            }
            var_by_dim <- numeric(length(dims_use))
            names(var_by_dim) <- rownames(loadings_mat[dims_use,])

            # determine features to use
            if (is.null(feats_use)) {
              feats_use <- colnames(input_data)[1:ncol(input_data)]
            } else if (is.integer(feats_use)) {
              feats_use <- colnames(input_data)[feats_use]
            }
            input_data <- input_data[,feats_use]
            # calculate total var
            mean_mat <- matrixStats::colMeans2(input_data)
            centered_mat <- input_data - mean_mat
            total_var = sum(centered_mat^2)

            # subset loadings mat
            loadings_mat <- loadings_mat[ , feats_use]

            # calculate variance explained per each individual dim
            for (d in dims_use) {
              # TODO: make separate function fo reconstructing data, separate for calculating r-squared
              reconst_d <- embeddings_mat[rownames(input_data), d, drop = FALSE] %*% loadings_mat[d, , drop = FALSE]
              resid_d <- input_data - reconst_d
              r2_d <- 1 - sum(resid_d^2)/total_var
              if (r2_d < 0) {
                # floor r2 at 0. while r2 can be negative, can't have lower than 0% of variance in data explained
                r2_d <- 0
              }
              var_by_dim[d] <- r2_d
            }
            # calculate total variance explained by dims picked
            reconst_all <- embeddings_mat[rownames(input_data), dims_use, drop = FALSE] %*% loadings_mat[dims_use, , drop = FALSE]
            resid_all <- input_data - reconst_all
            r2_all <- 1 - sum(resid_all^2)/total_var
            if (r2_all < 0) {
              r2_all <- 0
            }
            return(list(r2_by_dim = var_by_dim, r2_all = r2_all))
          })

# Methods for Running Dimensionality Reduction --------------------------------

#' Run PCA on assay data
#'
#' @param x = SummExpDR, SummarizedExperiment, or matrix.
#' matrix must have the form n features x m samples
#' @param i = assay to pull
#' @param suffix suffix to add, if SummExpDR. By default, reduced dim key is just 'PCA', but with suffix added it's 'PCA{suffix}'
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
            SampleID <- rownames(pca_coords)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            row_data <- S4Vectors::DataFrame(t(loadings_mat))
            # row_data$VarExplained <- pca_var/sum(pca_var)
            # recall that SummarizedExperiment objects have feats in rows and samples in cols
            pca_summ_exp <- create_FactorizedDR(assays = list(pca = t(pca_coords)), rowData = row_data, colData = col_data)
            return(pca_summ_exp)
          })


setMethod('runPCA',
          c('SummExpDR'),
          function(x, i, suffix = NULL, std_norm = TRUE) {
            if (is.null(suffix)) {
              suffix <- ''
            }
            if (is.numeric(i)) {
              assay_name <- SummarizedExperiment::assays(x@summ_exp)[i]
            } else {
              assay_name <- i
            }
            if (std_norm) {
              new_assay_name <- paste0(assay_name, '_std_norm')
              SummarizedExperiment::assay(x@summ_exp, new_assay_name) <- t(scale(t(SummarizedExperiment::assay(x@summ_exp, assay_name))))
              assay_use <- new_assay_name
            } else {
              assay_use <- assay_name
            }
            pca_summ_exp <- runPCA(SummarizedExperiment::assay(x@summ_exp, assay_use), std_norm = FALSE)
            # record assay used
            pca_summ_exp <- setSourceAssay(pca_summ_exp, assay_use)
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
#' @value returns SummExpDR object with
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
            tryCatch(stopifnot(is.character(assay_used)), error = function(e) {
              stop('assay not specified for Linear DR selected')
            })
            summ_exp <- getSummExp(x)
            assay_names <- SummarizedExperiment::assayNames(summ_exp)
            tryCatch(stopifnot(assay_used %in% assay_names), error = function(e) {
              stop('assay not in SummExpDR\'s SummarizedExperiment object\'s list of assay names')
            })
            # coords are in m samples x p features format
            # vmax loadings in p features x n embedding dims format
            # we transpose these in constructing new LinearDR (SummarizedExperiment) object
            data_coords <- t(SummarizedExperiment::assay(summ_exp, assay_used))
            new_coords <- data_coords %*% vmax_loadings
            rownames(new_coords) <- rownames(data_coords)
            row_data = S4Vectors::DataFrame(as.data.frame(t(vmax_loadings)))
            colnames(row_data) <- paste0(colnames(row_data), '_loading')
            SampleID <- rownames(new_coords)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            VmaxDR <- create_LinearDR(assays = list(varimax = t(new_coords)),
                                      rowData = row_data,
                                      colData = col_data)
            vmax_key <- paste0('varimax', suffix)
            x <- setReducedDims(x, key = vmax_key, value = VmaxDR)
            return(x)
          })


#' Run PCA on PCA (multistage PCA')
#'
#' Intended as a solution to capturing major sources of signal in multiple datasets while preventing
#' signal from being diminished due to presence of more features in data.
#' @include MultiExpIntegration.R
#' @param x multiExp object
#' @param assay_use assay to pull from
#' @param min_var = minimum pct variance explained per PC in PCA on each dataset
#' @param max_dims maximum number of dimensions to pull out
#' @value multiExp object with mPCA dim reduction added to reducedDim list
#' @export

setGeneric('mPCA', function(x, assay_use = 'imputed_mat', min_var = 0.05, max_dims = 5) standardGeneric('mPCA'))

setMethod('mPCA',
          signature = 'multiExp',
          function(x, assay_use = 'imputed_mat', min_var = 0.05, max_dims = 5) {
            # data is standard normalized
            x <- scaleExpData(x, assay_use, use_num_feats = FALSE)
            row_data <- rowData(x)
            unique_expts <- unique(row_data$expt)

            loadings_mats <- list()

            for (i in 1:length(unique_expts)) {
              expt <- unique_expts[i]
              print(paste('Running PCA for datatype', i, 'of', length(unique_expts)))
              expt_feats <- row_data[row_data$expt == expt , 'feat_id']
              # subset by features of given experiment, run PCA. data is already standard normalized
              x_subs <- subsetData(x, rows = expt_feats)
              x_subs <- runPCA(x_subs, i = 'scaled', std_norm = FALSE, suffix = '')
              var_expl_i <- varianceExplained(x_subs, 'PCA', dims_use = 1:max_dims)$r2_by_dim
              if (any(var_expl_i < min_var)) {
                if (all(var_expl_i < min_var)) {
                  warning('no PC met variance explained criterion, using first PC')
                  max_dims <- 1
                } else {
                  max_dims <- max(which(var_expl_i >= min_var))
                }
              }
              # loadings for datatype i
              V_i <- t(getLoadings(x_subs, 'PCA', paste0('PC', 1:max_dims)))
              # rescaling factors (1/sqrt(number of PCs)). do as a vector until
              # end of loop, then call diag to make diagonal matrix
              D_i <- rep(1/sqrt(ncol(V_i)), ncol(V_i))
              colnames(V_i) <- paste(expt, colnames(V_i), sep = '_')
              if (i == 1) {
                V <- V_i
                D <- D_i
              } else {
                missing_cols_x <- setdiff(colnames(V), colnames(V_i))
                missing_cols_y <- setdiff(colnames(V_i), colnames(V))
                pad_zeros <- function(mat, col_names) {
                  mat_append <- matrix(data = rep(0, nrow(mat)*ncol(mat)), nrow = nrow(mat), ncol = length(col_names))
                  rownames(mat_append) <- rownames(mat)
                  colnames(mat_append) <- col_names
                  mat <- cbind(mat, mat_append)
                }
                if (length(missing_cols_x) > 0) {
                  V_i <- pad_zeros(V_i, missing_cols_x)
                }
                if (length(missing_cols_y) > 0) {
                  V <- pad_zeros(V, missing_cols_y)
                }
                V <- rbind(V, V_i[,colnames(V)])
                D <- c(D, D_i)
                # gc(full = TRUE)
              }
            }
            D <- diag(D)
            # run PCA on PCA coords from PCA obtained in prior loop, scaled by # of PCs
            data_mat <- t(assay(x, 'scaled'))
            S <- V %*% D
            rot_coords1 <- data_mat %*% S
            P_res <- prcomp(S, scale. = FALSE, center = FALSE)
            Q <- P_res$rotation
            # final rotation matrix is multiplication of initial rotation matrix by scaling matrix by rotation matrix of PCA on PCA
            # final rotation matrix is rescaled to unit vectors. we do the unit normalization to simplify calculation of projection
            A <- S %*% Q
            A_len <- sqrt(t(matrixStats::colSums2(A^2)))
            A <- A/matrix(rep(A_len, nrow(A)), nrow = nrow(A), byrow = TRUE)
            colnames(A) <- paste0('mPC', 1:ncol(A))
            rownames(A) <- colnames(data_mat)
            # final coordinates
            P <- data_mat %*% A
            # set up return object
            SampleID <- rownames(P)
            col_data <- S4Vectors::DataFrame(SampleID = SampleID, row.names = SampleID)
            rownames(A) <- paste0(rownames(A), '_loading')
            row_data <- S4Vectors::DataFrame(t(A))
            mPCA_summ_exp <- create_FactorizedDR(assays = list(mPCA = t(P)), rowData = row_data, colData = col_data)
            mPCA_summ_exp <- setSourceAssay(mPCA_summ_exp, 'scaled')
            x <- setReducedDims(x, key = 'mPCA', value = mPCA_summ_exp)
            return(x)
          })
