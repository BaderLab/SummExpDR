###############################################################################
### Modified Legacy code from su2cproj package

###############################################################################
### Clustering Utilities, iterative clustering ###

### Overview: End goal is to iteratively run clustering (i_e_ run multiple times) with and without
###           resampling, assess clustering stability (with respect to
###           resampling), and for given value k (number of clusters)
###           report the clustering solution that best matches the solutions
###           found when resampling
###           Note that this general paradigm should work for any generic
###           clustering method, such as hierarchical clustering, spectral
###           clustering, kmeans, and apcluster_ For now, we'll just be using
###           spectral clustering for consistency with my previous work_
###           If we want to do NMF based clustering,
###           then the architecture of the clust_results slot
###           will have to be changed for this, since NMF package already has
###           a suite of tools for resampling based clustering on learned
###           factors, and reporting a final solution


#' Object to run iterative clustering/resampling procedure for single k
#'
#' @slot clusts_full NULL or matrix of clustering solutions obtained on
#' full dataset_
#' @slot clusts_rs NULL or matrix of clustering solutions obtained on
#' dataset when resampled in each iteration
#' @slot clusts_rs_perm NULL or matrix resulting from permuting rows of
#' clusts_rs
#' @slot sim_mat_full NULL or pairwise ARI matrix calculated by get_sim_mat on
#' clusts_full
#' @slot sim_mat_rs NULL or pairwise ARI_matrix calculated by get_sim_mat
#' on clusts_rs
#' @slot sim_mat_rs_perm NULL or pairwise ARI matrix calculated by
#' get_sim_mat on clusts_rs_perm
#' @slot unique_solns NULL or matrix containing unique solutions found in
#' clusts_rs by unique_solutions method
#' @slot soln_metadata NULL or data.frame containing metadata for unique_solns
#' solutions_ Metadata include proportion of solutions in clusts_full
#' accounted for by given unique solution, median similarity (ARI) to solutions
#' in clusts_rs, median absolute deviation (MAD) for similarity to solutions
#' in cluts_rs (derived from compar_mat_full_rs), and whether this cluster
#' solution is to be reported_ Currently, the decision to report is done first
#' by maximum median similarity to resampled solutions, and in case of ties,
#' by highest proportion of solutions accounted for in clusts_full_ In
#' case of further ties, pick the first solution out of the set of
#' solutions for which we need to break ties_
#' @slot compar_mat_full_rs NULL, or matrix M with M[i,j] = ARI
#' between solution i (from unique solutions) and solution j
#' (from clusts_rs)
#' @slot sim_rs_v_perm the results of a
#' wilcoxon ranksum test comparing distributions of pairwise ARI for
#' cluster labels obtained from resampled clustering and those obtained
#' from permuting resampled clustering based labels_
#' @export
#'
#' @examples
setClass('iter_clust_obj',
         slots = c(clusts_full = c('ANY'),
                   clusts_rs = c('ANY'),
                   clusts_rs_perm = c('ANY'),
                   sim_mat_full = c('ANY'),
                   sim_mat_rs = c('ANY'),
                   sim_mat_rs_perm = c('ANY'),
                   unique_solns = c('ANY'),
                   soln_metadata = c('ANY'),
                   compar_mat_full_rs = c('ANY'),
                   sim_rs_v_perm = c('ANY'))
)

## check validity of iter_clust_obj
check_iter_clust_obj <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  obj_slots <- slotNames(object)
  msg <- ''
  for (obj_slot_name in obj_slots){
    slot_val <- slot(object, obj_slot_name)
    if (obj_slot_name == 'soln_metadata') {
      if (is.null(slot_val)) {
        next
      } else if (!is.data.frame(slot_val)) {
        msg <- paste(msg,
                     'obj@soln_metadata must be dataframe',
                     'if not NULL',
                     sep = '\n')
      }
    } else if (obj_slot_name == 'sim_rs_v_perm') {
      if (is.null(slot_val)) {
        next
      } else if (!(is.list(slot_val))) {
        msg <- paste(msg,
                     'obj@sim_rs_v_perm must be list',
                     'if not NULL',
                     sep = '\n')
      }
    } else {
      if (is.null(slot_val) || is.matrix(slot_val)) {
        next
      } else {
        msg <- paste(msg,
                     paste0('obj@',obj_slot_name,' must be matrix'),
                     'if not NULL',
                     sep = '\n')
      }
    }
  }
  if (!msg == '') {
    stop(paste('Invalid clust_obj instance', msg, sep = '\n'))
  }
}

setValidity('iter_clust_obj', method = check_iter_clust_obj)

#' Constructor Function for iter_clust_obj
#'
#' @return new instance of iter_clust_obj
#' @export
#'
#' @examples
make_iter_clust_obj <- function() {
  new_obj <- new('iter_clust_obj',
                 clusts_full = NULL,
                 clusts_rs = NULL,
                 clusts_rs_perm = NULL,
                 sim_mat_full = NULL,
                 sim_mat_rs = NULL,
                 sim_mat_rs_perm = NULL,
                 unique_solns = NULL,
                 soln_metadata = NULL,
                 compar_mat_full_rs = NULL,
                 sim_rs_v_perm = NULL)
  return(new_obj)
}


#' make input matrix for clustering
#'
#' @param x matrix, with samples in rows and features in columns_
#' Note that in future we may try to incorporate SNF into this (i_e_
#' multiple data matrices input for SNF pipeline)
#' @param mode mode by which to construct affinity matrix_ if snf_affi,
#' then construct an affinity matrix using dist2 and affinityMatrix functions,
#' passing euclidean distance matrix to affinityMatrix function_
#' If std_norm is true, then data standard normalized prior to affinity matrix
#' construction_
#' If mxn, then return mxn matrix, standard normalized if std_norm is true_
#' Currently, only snf_affi and mxn supported
#' (i_e_ for hierarchical clustering)_
#' @param affi_K integer of length 1_ K parameter for \code{\link{SNFtool::affinityMatrix}}
#' @param sigma numeric of length 1_ sigma parameter for affinityMatrix
#' @param std_norm logical, length 1_ standard normalize data (normalize columns)?
#' @param as_kernMat logical, length 1_ coerce output matrix to class \code{\link{kernlab::kernelMatrix}}?
#' @return an affinity matrix
#' @examples
make_input_mat <- function(x, std_norm, mode = 'snf_affi', affi_K = 20, sigma = 0.5, as_kernMat = F) {
  if (nrow(x) == ncol(x)) {
    if (all(rownames(x) == colnames(x))) {
      warning('rownames and colnames of input x are equivalent_ Expecting
              mxn sample x feature matrix')
    }
  }
  if (std_norm) {
    x <- SNFtool::standardNormalization(x)
  }
  if (mode == 'snf_affi') {
    ## calculate Euclidean Distance

    euclid <- SNFtool::dist2(x, x)^(1/2)
    ## generate affinity matrix W as per Wang et al_ 2014
    inp_mat <- SNFtool::affinityMatrix(diff = euclid,
                                       K = affi_K,
                                       sigma = sigma)
  } else if (mode == 'mxn') {
    inp_mat <- x

  } else {
    stop(paste('mode', mode, 'not implemented'))
  }
  if (as_kernMat) {
    inp_mat <- kernlab::as.kernelMatrix(inp_mat)
  }
  return(inp_mat)
}

#' Run one iteration of clustering
#'
#' @param inp_mat matrix with numeric entries_ Affinity matrix if algorithm is spectral_kernlab or spectral_SNFtool,
#' sample x feature matrix normalized by columns if kmeans_ see \code{\link{make_inp_mat}}
#' @param algorithm character_ spectral_kernlab runs spectral clustering as implemented in kernlab package,
#' spectral_SNFtool runs spectral clustering as implemented in SNFtool package, kmeans runs
#' kmeans clustering as implemented in \code{\link{stats::kmeans}} package
#'
#' @return numeric_ cluster labels in the order of samples in the input matrix
#' @export
#'
#' @examples
clust_once <- function(inp_mat, algorithm, k, seed) {
  set.seed(seed)
  if (algorithm == 'spectral_SNFtool') {
    clust_labs <- SNFtool::spectralClustering(inp_mat,
                                              K = k,
                                              type = 3)
  } else if (algorithm == 'spectral_kernlab') {
    clust_labs <- kernlab::specc(x = inp_mat, centers = k)@.Data
  } else if (algorithm == 'kmeans') {
    clust_labs <- stats::kmeans(inp_mat, centers = k)
    clust_labs <- clust_labs$cluster
  } else {
    stop(paste('algorithm', algorithm, 'not yet included for function'))
  }
  return(clust_labs)
}

#' Iteratively run clustering w/ or w/o resampling
#' @inheritParams clust_once
#' @param x matrix_ samples in rows, features in columns
#' in future versions may allow for a list of matrices
#' for SNF based clustering
#' @param num_iter numeric, length 1_ number of iterations of clustering_
#' Note that if we are running a deterministic clustering algorithm
#' without resampling, this argument automatically gets set to 1
#' @param algorithm
#' @param resample logical, length 1_ whether or not to resample data
#' in each iteration_ If this is the case, and an affinity matrix must
#' be calculated (this is the case for spectral clustering), then
#' the affinity matrix is recalculated for each iteration based on
#' data in resample_
#' @param pct_resample percentage of data to resample_ Specify as number
#' greater than 0 and less than 1
#' @param affi_K numeric, length 1_ K value for affinityMatrix_
#' used if algorithm = 'spectral_SNFtool' or 'spectral_kernlab'
#' @param sigma numeric of length 1_ sigma value for affinityMatrix_
#' used if algorithm = 'spectral_kernlab' or 'spectral_SNFtool'
#' @param k integer of length 1_ number of clusters
#' @param base_seed integer of length 1_ lowest seed in range
#' base_seed:(base_seed + num_iter - 1) used for resampling,
#' and in case of clustering with restarts, clustering_ For resampling,
#' in each iteration, the seed used is incremented up in this range,
#' and before clustering the seed is reset to the current seed in the range
#' (i_e_ in a given iteration the seed set for resampling and clustering are
#' the same)
#' @param std_norm standard normalize data prior to clustering?
#' and is set prior to sampling the data_
#' @param n_cores number of cores to use
#' @param pbar progress bar
#' @return a matrix where each row is a clustering result_ If resampling
#' applied, then each clustering result returned in the matrix has cluster
#' labels for samples within the resampling subset of the given iteration,
#' and zero labels for samples not included in the resample_
#' @details TODO
#' @export
#'
#' @examples
iter_clust <- function(x,
                       k,
                       std_norm = T,
                       num_iter = 200,
                       algorithm = 'spectral_kernlab',
                       resample = F,
                       pct_resample = 0.8,
                       affi_K = 20,
                       sigma = 0.5,
                       base_seed = 2,
                       n_cores = NULL,
                       pbar = FALSE) {
  # print('running iter_clust')
  num_samples <- nrow(x)
  num_resample <- floor(num_samples*pct_resample)
  sample_names <- rownames(x)
  if (algorithm %in% c('spectral_SNFtool', 'spectral_kernlab')) {
    inp_mat_mode <- 'snf_affi'
  } else if (algorithm == 'kmeans') {
    inp_mat_mode <- 'mxn'
  } else {
    stop(paste('algorithm', algorithm, 'not included for this function'))
  }

  ## if input matirx must be of class kernelMatrix, set it to be one
  ## in make_input_mat
  if (algorithm == 'spectral_kernlab') {
    as_kernMat <- T
  } else {
    as_kernMat <- F
  }
  if (!resample) {
    # print('not resampling, using predetermined input matrix')
    inp_mat <- make_input_mat(x,
                              std_norm = std_norm,
                              mode = inp_mat_mode,
                              affi_K = affi_K,
                              sigma = sigma,
                              as_kernMat = as_kernMat)
  } else {
    # print('resampling')
  }
  output_mat <- matrix(nrow = 0, ncol = num_samples)
  # colnames(output_mat)
  ## if clustering algorithm used is deterministic (i_e_ does not change with different seeds),
  ## set iterations to 1 if not resampling
  if (algorithm == 'spectral_SNFtool') {
    if (!resample) {
      num_iter <- 1
    }
  }
  do_iter <- function(i) {
    seed_use <- base_seed + i - 1
    if (resample) {
      set.seed(seed_use)
      sample_inds <- sample(1:num_samples, num_resample, replace = F)
      # subset of samples to use
      samples_use <- sample_names[sample_inds]
      inp_mat <- make_input_mat(matrix(x[samples_use,], nrow = length(samples_use), ncol = ncol(x)),
                                std_norm = std_norm,
                                mode = inp_mat_mode,
                                affi_K = affi_K,
                                sigma = sigma,
                                as_kernMat = as_kernMat)
      clust_result_rs <- clust_once(inp_mat,
                                    algorithm = algorithm,
                                    k = k,
                                    seed = seed_use)
      # make new vector with all samples, with cluster labels for missing samples in resampling set to 0
      names(clust_result_rs) <- samples_use
      clust_result <- numeric(num_samples)
      names(clust_result) <- sample_names
      clust_result[samples_use] <- clust_result_rs
    } else {
      clust_result <- clust_once(inp_mat,
                                 algorithm = algorithm,
                                 k = k,
                                 seed = seed_use)
      # names(clust_result) <- sample_names
    }
    result_mat <- matrix(clust_result, nrow = 1)
    colnames(result_mat) <- sample_names
    return(as.data.frame(result_mat))
  }
  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores()
    message(paste('automatically assigning', n_cores, 'cores for clustering'))
  }
  if (n_cores > num_iter) {
    warning(paste('reducing number of cores to num_iter as n_cores > num_iter'))
    n_cores <- num_iter
  }
  if (n_cores == 1) {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar, exportglobals = TRUE)
  }
  output_list <- BiocParallel::bplapply(1:num_iter, FUN = function(x) {do_iter(x)}, BPPARAM = bpparam)
  output_mat <- as.matrix(dplyr::bind_rows(output_list))
  # for (i in 1:num_iter) {
  #   seed_use <- base_seed + i - 1
  #   if (resample) {
  #     set.seed(seed_use)
  #     sample_inds <- sample(1:num_samples, num_resample, replace = F)
  #     samples_use <- sample_names[sample_inds]
  #     inp_mat <- make_input_mat(matrix(x[samples_use,], nrow = length(samples_use), ncol = ncol(x)),
  #                               std_norm = std_norm,
  #                               mode = inp_mat_mode,
  #                               affi_K = affi_K,
  #                               sigma = sigma,
  #                               as_kernMat = as_kernMat)
  #     clust_result_rs <- clust_once(inp_mat,
  #                                   algorithm = algorithm,
  #                                   k = k,
  #                                   seed = seed_use)
  #     names(clust_result_rs) <- samples_use
  #     clust_result <- numeric(num_samples)
  #     names(clust_result) <- sample_names
  #     clust_result[samples_use] <- clust_result_rs
  #   } else {
  #     clust_result <- clust_once(inp_mat,
  #                                algorithm = algorithm,
  #                                k = k,
  #                                seed = seed_use)
  #     names(clust_result) <- sample_names
  #   }
  #   output_mat <- rbind(output_mat, clust_result)
  # }
  ## name clustering solutions by iteration
  rownames(output_mat) <- 1:num_iter
  return(output_mat)
}
#' Run iterative clustering w/ and w/o resampling, permutation
#' @description run multiple rounds of clustering on data in
#' iter_clust_obj, doing so on full dataset without resampling,
#' and again on dataset with resampling applied, and obtain
#' permuted labels for cluster labels obtained from resampling
#' @inheritParams  iter_clust
#' @param object
#' @param inp_data
#' @param k
#' @param num_iter
#' @param std_norm
#' @param algorithm
#' @param pct_resample
#' @param affi_K
#' @param sigma
#' @param base_seed
#' @param ... other args for iter_clust
#' @details iter_clust is run first on the whole dataset (i_e_
#' no resampling), generating a matrix with num_iter clustering solutions
#' with labels for all samples_ This result goes to the clusts_full slot
#' of the input object_ iter_clust is then run again, this time with
#' resampling in each iteration, generating a similar matrix of clustering
#' solutions, but for each solution labels are only given to samples that
#' were included in the sampling for the given iteration_ This result is
#' put in the slot clusts_rs_ We then generate a copy of the matrix
#' given to clusts_rs and permute the labels for each solution, effectively
#' given num_iter 'solutions' with random labels_ The permuted labels matrix
#' is assigned to the slot clusts_rs_perm
#' @return iter_clust_obj, with slots clusts_full, clusts_rs,
#' and clusts_rs_perm determined_
#' @export
#'
#' @examples

do_clust <- function(object,
                     inp_data,
                     k,
                     num_iter = 200,
                     std_norm = T,
                     algorithm = 'spectral_kernlab',
                     pct_resample = 0.8,
                     affi_K = 20,
                     sigma = 0.5,
                     base_seed = 2,
                     ...) {
  clusts_full <- iter_clust(x = inp_data,
                            k = k,
                            num_iter = num_iter,
                            std_norm = std_norm,
                            algorithm = algorithm,
                            resample = F,
                            pct_resample = pct_resample,
                            affi_K = affi_K,
                            sigma = sigma,
                            base_seed = base_seed)
  clusts_rs <- iter_clust(x = inp_data,
                          k = k,
                          num_iter = num_iter,
                          std_norm = std_norm,
                          algorithm = algorithm,
                          resample = T,
                          pct_resample = pct_resample,
                          affi_K = affi_K,
                          sigma = sigma,
                          base_seed = base_seed,
                          ...)
  clusts_rs_perm <- perm_labels(M = clusts_rs,
                                base_seed = base_seed)
  object@clusts_full <- clusts_full
  object@clusts_rs <- clusts_rs
  object@clusts_rs_perm <- clusts_rs_perm
  check_iter_clust_obj(object)
  return(object)
}


#' Compute pairwise adjusted rand index matrix
#'
#' @param M matrix of clustering solutions, each row being a different
#' solution_ must have > 1 row
#' @description compute pairwise adjusted rand index
#' (ARI) between all clustering solutions in input matrix
#' @details
#' @return sim_mat matrix of pairwise similarities for input clustering
#' solutions_ Let sim_mat be the output matrix_ then sim_mat[i,j] is the ARI
#' between solutions at rows i and j of input matrix x_ If solutions i and j
#' do not have a completely overlapping set of samples (due to resampling
#' applied in iter_clust), then ARI is calculated between all samples in the
#' intersection of samples between the two solutions
#' @export
#'
#' @examples
get_sim_mat <- function(M) {
  iterations <- nrow(M)
  if (iterations == 1) {
    stop('nrow(M) (clustering solutions) must be > 1')
  } else {
    sim_mat <- matrix(rep(0, iterations^2),
                      nrow = iterations,
                      ncol = iterations)
    ## do all pairs s_t_ i != j
    for (i in 1:(iterations - 1)) {
      for (j in (i+1):iterations) {
        Mi <- M[i,]
        Mj <- M[j,]
        if (any(Mi == 0) || any(Mj == 0)) {
          ## if we are dealing with zero indices (i_e_ indices without labels),
          ## only compare samples that are common between the two clustering
          ## solutions
          zero_ind_i <- which(Mi == 0)
          zero_ind_j <- which(Mj == 0)
          inds_remove <- union(zero_ind_i, zero_ind_j)
          Mi <- Mi[-inds_remove]
          Mj <- Mj[-inds_remove]
        }
        sim_mat[i,j] <- mclust::adjustedRandIndex(Mi, Mj)
      }
    }
    ## make matrix symmetric, add diagonal since ARI of set of labels
    ## with self is 1
    sim_mat <- sim_mat + t(sim_mat)
    sim_mat <- sim_mat + diag(x = 1, nrow = iterations, ncol = iterations)
  }
  return(sim_mat)
}

#' permute cluster labels in matrix
#'
#' @param M matrix output by iter_clust
#' @param base_seed smallest seed used_ seeds used
#' for permutations range from base_seed to base_seed
#' plus nrow(M) (i_e_ number of iterations of clustering)
#' minus 1
#' @return matrix, with each row permuted
#' @export
#'
#' @examples
perm_labels <- function(M, base_seed = 2)
{
  for (i in 1:nrow(M)) {
    set.seed((base_seed + i - 1)^2)
    M[i,] <- gtools::permute(M[i,])
  }
  return(M)
}

#' Calculate pairwise ARI matrices
#' @description Calculate pairwise ARI matrices for resampled clustering
#' solutions, permuted labels, and, if multiple funs for clustering
#' on full dataset, for solutions found for full dataset_
#' @param object iter_clust_obj output by do_clust_
#'
#' @return
#' @export
#'
#' @examples
calc_sim_mats <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  if (any(c(is.null(object@clusts_full),
            is.null(object@clusts_rs),
            is.null(object@clusts_rs_perm)
  )
  )
  ) {
    stop('null values for 1 or more of clustering slots_ likely that
         do_clust was not run_')
  }
  if (nrow(object@clusts_full) == 1) {
    object@sim_mat_full <- NULL
  } else if (nrow(object@clusts_full) > 1) {
    object@sim_mat_full <- get_sim_mat(object@clusts_full)
  } else {
    stop('expect more than one row for obj@clusts_full')
  }
  object@sim_mat_rs <- get_sim_mat(object@clusts_rs)
  object@sim_mat_rs_perm <- get_sim_mat(object@clusts_rs_perm)
  check_iter_clust_obj(object)
  return(object)
}

#' return vector of pairwise similarities
#'
#' @param sim_mat output of function get_sim_mat
#' @details returns all pairwise ARI from sim_mat,
#' i_e_ sim_mat[i,j] for all i, j, j > i (avoid counting twice)
#' @return numeric vector of all non-self pairwise ARI
#' @export
#'
#' @examples
get_sim_dist <- function(sim_mat) {
  vect_return <- numeric(0)
  for (i in 1:(nrow(sim_mat) - 1)) {
    for (j in (i+1):ncol(sim_mat)) {
      vect_return <- c(vect_return, sim_mat[i,j])
    }
  }
  return(vect_return)
}

#' Find unique solutions in full dataset clustering
#' and proportion of solutions accounted for
#'
#' @param object iter_clust_obj output by calc_sim_mats
#'
#' @return iter_clust_obj with unique solutions specified and
#' soln_metadata slot initialized as a data.frame containing
#' the proportion of solutions in the full dataset clusterings
#' accounted for by each unique solution_
#' @export
#'
#' @examples
unique_solutions <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  clusts_full <- object@clusts_full
  if (is.null(clusts_full)) {
    stop('clusts_full set to NULL_ likely that do_clust was not run
         for this object')
  }
  ## make empty data.frame to contain metadata for solutions
  soln_metadata <- data.frame(solution = numeric(0),
                              proportion = numeric(0))

  ## find unique solutions
  unique_solns <- matrix(nrow = 0, ncol = ncol(clusts_full))
  colnames(unique_solns) <- colnames(clusts_full)
  if (nrow(clusts_full) == 1) {
    ## we only ran 1 round of clustering, so therefore 1 solution
    clust_soln_i <- as.numeric(clusts_full[1,])
    unique_solns <- rbind(unique_solns, clust_soln_i)
    new_row <- data.frame(solution = as.integer(1), proportion = 1)
    soln_metadata <- rbind(soln_metadata, new_row)
  } else if (nrow(clusts_full > 1)) {
    ## more than 1 round of clustering (e_g kmeans with restarts)
    sim_mat_full <- object@sim_mat_full
    if (is.null(sim_mat_full)) {
      stop('more than 1 round of clustering performed, need to determine
           similarity of solutions (use calc_sim_mats(object)) prior to
           determining unique solutions')
    }
    ## numeric vector to store indices of solutions for which
    ## an identical solution has already been found
    found_inds <- integer(0)
    ## records solution number
    soln_num <- as.integer(1)
    ## procedure is to iterate through rows, check if an identical
    ## solution has already been found (by checking found_inds),
    ## and if not, find the solution corresponding to the current
    ## row (row corresponds to iteration in clusts_full)
    for (i in 1:nrow(clusts_full)) {
      if (!(i %in% found_inds)) {
        clust_soln_i <- clusts_full[i,]
        unique_solns <- rbind(unique_solns, clust_soln_i)
        is_equiv_soln <- sim_mat_full[i,] == 1
        found_inds <- c(found_inds, which(is_equiv_soln))
        prop_soln <- sum(as.numeric(is_equiv_soln))/ncol(sim_mat_full)
        new_row <- data.frame(solution = soln_num, proportion = prop_soln)
        soln_metadata <- rbind(soln_metadata, new_row)
        soln_num <- as.integer(soln_num + 1)
      } else {
        next
      }
    }
  } else {
    stop('nrow(clusts_full_rs must be > 0')
  }
  colnames(unique_solns) <- colnames(clusts_full)
  object@unique_solns <- unique_solns
  object@soln_metadata <- soln_metadata
  check_iter_clust_obj(object)
  return(object)
}

#' Compare distributions of pairwise ARI for resampled cluster solutions
#' and permuted labels
#'
#' @param object iter_clust_obj
#'
#' @return iter_clust_obj with sim_rs_v_perm filled in with
#' result of one sided wilcoxon ranksum test_ Note that
#' the alternative hypothesis is that the distribution of
#' ARI values for the resampled clustering solutions is
#' higher than that of the the permuted labels
#' @export
#'
#' @examples
compare_rs_perm <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  if (is.null(object@sim_mat_rs) || is.null(object@sim_mat_rs_perm)) {
    stop('Pairwise ARI for resampled clustering solutions
         and permuted labels must be calculated_ use
         calc_sim_mats(iter_clust_obj) to do this')
  }
  rs_ARI_vals <- get_sim_dist(object@sim_mat_rs)
  rs_perm_ARI_vals <- get_sim_dist(object@sim_mat_rs_perm)
  wcox_ranksum <- wilcox.test(x = rs_ARI_vals,
                              y = rs_perm_ARI_vals,
                              alternative = 'greater')
  object@sim_rs_v_perm <- wcox_ranksum
  check_iter_clust_obj(object)
  return(object)
}

#' Compare clustering solutions on full data to clustering
#' solutions on resampled data to determine similarity
#'
#' @param object iter_clust_obj
#'
#' @return iter_clust_obj with compar_mat_full_rs assigned as a matrix
#' of pairwise similarities between all clustering solutions on full
#' dataset (rows) with all resampled clustering solutions, and with
#' soln_metadata slot having a column added for median ARI compared to
#' resampled solutions
#' @export
#'
#' @examples
compare_full_rs <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  if (is.null(object@clusts_full) || is.null(object@clusts_rs)) {
    stop('clusts_full set to NULL, clusts_rs set to NULL_
         Run do_clust(iter_clust_object) prior to running this
         function')
  } else if (is.null(object@unique_solns)) {
    stop('need to find unique solutions_ run
         iter_clust_obj <- calc_sim_mats(iter_clust_obj), then
         iter_clust_obj <- unique_solutions(iter_clust_obj),
         prior to running this function')
  }
  unique_solns <- object@unique_solns
  clusts_rs <- object@clusts_rs
  if (!all(colnames(unique_solns) == colnames(clusts_rs))) {
    stop('Expect that unique_solns and clusts_rs have identical
         composition and ordering of samples')
  }
  n_unique <- nrow(unique_solns)
  n_rs <- nrow(clusts_rs)
  compar_mat <- matrix(rep(NA, n_unique*n_rs),
                       nrow = n_unique,
                       ncol = n_rs)
  for (i in 1:n_unique) {
    for (j in 1:n_rs) {
      rs_soln_j <- clusts_rs[j,]
      inds_keep <- which(rs_soln_j != 0)
      unique_solns_i <- unique_solns[i,inds_keep]
      rs_soln_j <- rs_soln_j[inds_keep]
      compar_mat[i,j] <- mclust::adjustedRandIndex(unique_solns_i,
                                                   rs_soln_j)
    }
  }
  if (any(is.na(compar_mat))) {
    stop('expect all elements of comparison matrix to be of numeric
         value')
  }
  object@compar_mat_full_rs <- compar_mat
  check_iter_clust_obj(object)
  return(object)
}
#' Report solution to use
#'
#' @param object iter_clust_obj
#'
#' @return iter_clust_obj with soln_metadata slot having a column added
#' for median ARI compared to resampled solutions, and a column added
#' denoting whether or not the solution is the one to be reported
#' for this iter_clust_obj_ Note that in case of algorithms that
#' only have 1 clustering performed on full data (e_g_ SNFtool
#' implementation of spectral clustering), this degenerates
#' to picking the only solution_
#' @details if more than one unique solution to pick from,
#' pick the solution with the highest median ARI when compared to
#' resampled solutions_ In case of ties, pick solution with highest
#' proportion of clustering solutions in full dataset_ In cases
#' of further ties, pick the first solution our of the set of
#' the ties for proportion
#' @export
#'
#' @examples
pick_soln <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  compar_mat <- object@compar_mat_full_rs
  soln_metadata <- object@soln_metadata
  compar_mat <- object@compar_mat_full_rs
  if (is.null(soln_metadata)) {
    stop('expect that metadata for unique solutions exists_
         run iter_clust_obj <- unique_solutions(iter_clust_obj),
         then iter_clust_obj <- run compar_full_rs(iter_clust_obj)
         to get comparisons of unique solutions to solutions from
         resampled data')
  }
  else if (is.null(compar_mat)) {
    stop('expect that comparison of unique solutions to
         resampled solutions already done_ run compare_full_rs,
         prior to running this function')
  }
  ## get median ARI with resampled data for each unique solution
  med_ARI_rs <- apply(compar_mat, MARGIN = 1, FUN = median)
  ## pick solution to report
  soln_report <- logical(nrow(soln_metadata))
  max_med_ARI <- which(med_ARI_rs == max(med_ARI_rs))
  if (length(max_med_ARI) > 1) {
    ## break ties with proportion
    tied_solns <- max_med_ARI
    proportion <- soln_metadata$proportion
    max_prop_tied <- max(proportion[tied_solns])
    tie_winner <- intersect(tied_solns, which(proportion == max_prop_tied))
    if (length(tie_winner) > 1) {
      ## if more than one tiebreaking winner, simply use first of tie winners
      tie_winner <- tie_winner[1]
    }
    soln_report[tie_winner] <- T
  } else if (length(max_med_ARI) == 1) {
    soln_report[max_med_ARI] <- T
  } else {
    stop('should have at least one solution with max median ARI
         comparing to resampled solutions')
  }
  if (!sum(as.numeric(soln_report) == 1)) {
    stop('should have 1 and only 1 solution to report_ function has a bug')
  }
  soln_metadata$med_ARI_rs <- med_ARI_rs
  soln_metadata$soln_report <- soln_report
  object@soln_metadata <- soln_metadata
  check_iter_clust_obj(object)
  return(object)
}

#' Run do_clust, calc_sim_mats, unique_solutions, compare_rs_perm,
#' compare_full_rs, pick_soln in sequence
#'
#' @param inp_data mxn matrix of data
#' @inheritParams do_clust
#' @param ... other args for iter_clust
#' @return iter_clust_obj with a solution reported as well as associated
#' information
#' @details if more than one unique solution to pick from,
#' pick the solution with the highest median ARI when compared to
#' resampled solutions_ In case of ties, pick solution with highest
#' proportion of clustering solutions in full dataset_ In cases
#' of further ties, pick the first solution our of the set of
#' the ties for proportion
#' @export
#'
#' @examples
run_iter_clust_functions <- function( inp_data,
                                      k,
                                      num_iter = 200,
                                      std_norm = T,
                                      algorithm = 'spectral_kernlab',
                                      pct_resample = 0.80,
                                      affi_K = 20,
                                      sigma = 0.5,
                                      base_seed = 2,
                                      ...) {
  object <- make_iter_clust_obj()
  object <- do_clust(inp_data = inp_data,
                     object = object,
                     k = k,
                     num_iter = num_iter,
                     std_norm = std_norm,
                     algorithm = algorithm,
                     pct_resample = pct_resample,
                     affi_K = affi_K,
                     sigma = sigma,
                     base_seed = base_seed,
                     ...)
  object <- calc_sim_mats(object)
  object <- unique_solutions(object)
  object <- compare_rs_perm(object)
  object <- compare_full_rs(object)
  object <- pick_soln(object)
  check_iter_clust_obj(object)
  return(object)
}

#' Retrieve Metadata for Unique Solutions
#' @param object iter_clust_obj
#' @return data.frame, metadata for unique solutions_
#' For more info, see documentation for class iter_clust_obj
get_soln_metadata <- function(object) {
  stopifnot(is(object, 'iter_clust_obj'))
  return(object@soln_metadata)
}

#' Retrieve Solution From iter_clust_obj
#'
#' @param object iter_clust_object
#' @param  get_best logical of length 1_ if true, retrieve solution
#' automatically picked to be reported_ if flase, must set soln_num
#' to an integer value
#' @param soln_num numeric of length 1_ which solution to pick, if
#' get_best set to false_ must be set to an integer value within 1:
#' number of solutions if get_best set to false
#' @return named numeric vector of cluster labels for solution
#' selected
get_soln_iter_clust <- function(object, get_best = T, soln_num = NULL) {
  stopifnot(is(object, 'iter_clust_obj'))
  soln_metadata <- get_soln_metadata(object)
  unique_solns <- object@unique_solns
  if (get_best) {
    soln_use <- which(soln_metadata$soln_report == T)
    if (length(soln_use) != 1) {
      stop('bad solution metadata_ expect 1 and only 1
           solution picked to report')
    }
    soln_return <- unique_solns[soln_use,]
    names(soln_return) <- colnames(unique_solns)
  } else if (is.numeric(soln_num) && length(soln_num) == 1) {
    if (!soln_num %in% 1:nrow(unique_solns)) {
      stop(paste0('bad argument soln_num = ', soln_num, '_ Must be
                  betwen 1 and the number of unique solutions found,
                  (', nrow(unique_solns), ')'))
    }
    soln_return <- unique_solns[soln_num,]
    names(soln_return) <- colnames(unique_solns)
  } else {
    stop('get_best must be true or soln_num must be a numeric of length 1')
  }
  return(soln_return)
}

#' Report metrics for qualitiy control for single k
#' @param iter_clust_obj object of class iter_clust_obj, produced from
#' \code{\link{run_iter_clust_functions}}
#' @param thresh_pval if wilcoxon ranksum test comparing ARI for resampled clusterings
#' to that for permuted labels gives p_value <= this argument, consider
#' the result significant
#' @param thresh_ARI threshold ARI for pairwise comparisons of resampling based
#' clustering solutions, and pairwise comparisons of reported clustering solution
#' to resampling based clustering solutions
#' run_iter_clust_functions
#' @param thresh_proportion threshold proportion of
#' @details the following QC metrics will be reported:
#' statistical significance for wilcoxon rank sum test comparing pairwise ARI
#' for resampled clustering solutions and those obtained for permuted labels,
#' determined with specified p_value threshold; p_value for said test (is_signif);
#' proportion of pairwise ARI for resampling based clustering solutions with
#' ARI above specified threshold (prop_rs_v_rs);
#' proportion of pairwise ARI in comparisons of reported solution (for given k)
#' to resampling based clusters (prop_rpt_v_rs);
#' does proportion of pairwise ARI values for the reampled clustering comparisons
#' and, and that calculated for the reported solution x resampled solutions comparisons,
#' exceed a specified threshold proportion (meets_qc);
#' These values are all stored in teh data.frame that is returned
#' @return data.frame with 1 row, and the following columns:
#' - is_signif: logical
#' - prop_rs_v_rs: numeric
#' - prop_rpt_v_rs: numeric
#' - meets_qc: logical
#' see details for more info
#' @export

get_qc_k <- function(iter_clust_obj,
                     thresh_pval = 0.05,
                     thresh_ARI = 0.90,
                     thresh_proportion = 0.75) {
  # p-value from wilcoxon test
  p_val_pass <- iter_clust_obj@sim_rs_v_perm$p.value < thresh_pval

  sim_dist_rs_v_rs <- get_sim_dist(iter_clust_obj@sim_mat_rs)
  prop_rs_v_rs <- sum(sim_dist_rs_v_rs > thresh_ARI)/length(sim_dist_rs_v_rs)
  ARI_pass_rs_v_rs <- (prop_rs_v_rs > thresh_proportion)

  sim_dist_rpt_v_rs <- as.numeric(iter_clust_obj@compar_mat_full_rs)
  prop_rpt_v_rs <- sum(sim_dist_rpt_v_rs > thresh_ARI)/length(sim_dist_rpt_v_rs)
  ARI_pass_rpt_v_rs <- (prop_rpt_v_rs > thresh_proportion)

  meets_qc <- all(c(p_val_pass, ARI_pass_rs_v_rs, ARI_pass_rpt_v_rs))

  return(data.frame(is_signif = p_val_pass,
                    prop_rs_v_rs = prop_rs_v_rs,
                    prop_rpt_v_rs = prop_rpt_v_rs,
                    meets_qc = meets_qc))
}

###############################################################################
### Clustering Utilities, iterative clustering, multiple values k # clusters

### OVERIVEW: The functions presented in clustering_utils provide a framework
###           for obtaining cluster labels and assessing stability at 1 value
###           k (k = # of clusters). The set of functions presented here
###           intend to automate running the iter_clust_obj related funcions
###           for multiple values k, and picking the clustering solution
###           to report in a somewhat principled manner.


#' object to run clustering analyses for multiple values k
#' @slot single_k_analyses either NULL or a named list of iter_clust_obj
#' objects. Naming convetion for list: single_k_analyses[['k']] =
#' iter_clust_obj run for k clusters
#' @slot report data.frame with rows corresponding to clustering analyses
#' at k = (row #), and with first (n -1) columns corresponding to variables
#' of interest used to determine the value k selected.
#' Columns include:
#' - is.significant = are pairwise ARI for resampled clusters significantly
#' higher than would be expected under a null hypothesis of permuted labels?
#' - rs.ARI.gt.thresh = proportion of pairwise ARI greater than threshold set
#' - full.v.rs.ARI.gt.thresh = proportion of pairwise ARI between reported
#' solution on full data and resampled solutions exceeding threshold
#' -report.k = report solution for this value k?
setClass('multi_k_clust',
         slots = c('single_k_analyses' = c('ANY'),
                   'report' = c('ANY')))

valid_multi_k_clust <- function(object) {
  if (!is.null(object@single_k_analyses)) {
    single_k_analyses <- object@single_k_analyses
    if (!is.list(single_k_analyses)) {
      stop('object@single_k_analyses must be list')
    }
    for (i in 1:length(single_k_analyses)) {
      if (!is(single_k_analyses[[i]], 'iter_clust_obj')) {
        stop('elements of single_k_analyses must be of class iter_clust_obj')
      }
    }
  }
  if (!is.null(object@report)) {
    if (!is.data.frame(object@report)) {
      stop('report must be NULL or data.frame')
    }
  }
}

setValidity('multi_k_clust', valid_multi_k_clust)

#' create a multi_k_clust_object
#' @return \code{\link{multi_k_clust_object}}

create_multi_k_clust <- function() {
  new_obj <- new('multi_k_clust',
                 single_k_analyses = NULL,
                 report = NULL)
  return(new_obj)
}

#' run clustering/evaluation workflow for multiple k
#' @inheritParams do_clust
#' @param object multi_k_clust object
#' @param k_use integer vector, length >= 1. values k to perform
#' clustering for
#' @return multi_k_clust object with single.k.analysis slot assigned
#' as a list of iter_clust_obj objects for which clustering and
#' evaluation of full data and resampled solutions has been done.
#' Each element of the list is named according to the value k that
#' was specified for the numnber of clusters.
#' for more dtails, see run_iter_clust_functions.

run_clust_multi_k <- function(object,
                              inp_data = NULL,
                              k_use = 2:10,
                              num_iter = 200,
                              std_norm = NULL,
                              algorithm = 'spectral_kernlab',
                              pct_resample = 0.80,
                              affi_K = 20,
                              sigma = 0.5,
                              base_seed = 2,
                              ...) {
  if (!is.integer(k_use)) {
    stop('argument k_use must be integer')
  } else if (any(k_use == 1)) {
    stop('k_use must be specified as values k greater than 1')
  }
  if (!is.null(object@single_k_analyses)) {
    stop('object@single_k_analyses is not NULL. object may be output
         of this function')
  }
  single_k_analyses <- list()
  for (i in k_use) {
    iter_clust_i <- run_iter_clust_functions(inp_data = inp_data,
                                             num_iter = num_iter,
                                             std_norm = std_norm,
                                             algorithm = algorithm,
                                             k = i,
                                             pct_resample = pct_resample,
                                             affi_K = affi_K,
                                             sigma = sigma,
                                             base_seed = base_seed,
                                             ...)
    single_k_analyses[[as.character(i)]] <- iter_clust_i
  }
  object@single_k_analyses <- single_k_analyses
  valid_multi_k_clust(object)
  return(object)
}

#' Take Input Data through clustering with multiple values k
#'
#' @inheritParams do_clust
#' @param x either an mxn sample x feature matrix, or a SummExpDR object containing input data
#' or pca results
#' @param assay_use if SummExpDR object, which assay to use from expt.data slot
#' @param analysis_name name for clustering analysis. argument unused if x is a matrix
#' @param redDimKey if NULL, use coordinates
#' @param dims_use
#' @param std_norm standard normalize data prior to clustering. for spectral clustering, this
#' means standard normalization prior to construciton of affinity matrix. If unspecified,
#' will exhibit following default behavior: For multi_clust_k_obj, defaults to TRUE;
#' for SummExpDR, if redDimKey is NULL, defaults to TRUE, if redDimKey is non-NULL,
#' defaults to FALSE
#' @param k_use values k to use for clustering
#' @param num_iter
#' @param algorithm
#' @param pct_resample
#' @param affi_K
#' @param sigma
#' @param base_seed
#' @value if matrix input, \code{\link{multi_k_clust}} object with clustering run for all
#' specified values k in k_use, as well as qc metric determined. If a SummExpr object, a SummExpr object
#' with a multi_k_clust_obj added as entry 'clust.analysis' to analyses slot
#' @export

run_multi_k_functions  <- function(x,
                                   assay_use = 1,
                                   analysis_name = 'clust.analysis',
                                   redDimKey = F,
                                   dims_use = NULL,
                                   std_norm = NULL,
                                   k_use = 2:10,
                                   num_iter = 200,
                                   algorithm = 'spectral_kernlab',
                                   pct_resample = 0.80,
                                   affi_K = 20,
                                   sigma = 0.5,
                                   base_seed = 2,
                                   thresh_pval = 0.05,
                                   thresh_ARI = 0.90,
                                   thresh_proportion = 0.90,
                                   ...) {
  if (is(x, 'matrix')) {
    if (is.null(std_norm)) {
      std_norm <- T
    }
    multi.k.obj <- create_multi_k_clust()
    multi.k.obj <- run_clust_multi_k(object = multi.k.obj,
                                     inp_data = x,
                                     k_use = k_use,
                                     std_norm = std_norm,
                                     num_iter = num_iter,
                                     algorithm = algorithm,
                                     pct_resample = pct_resample,
                                     affi_K = affi_K,
                                     sigma = sigma,
                                     base_seed = base_seed,
                                     ...)
    multi.k.obj <- get_qc_all_k(multi_k_clust = multi.k.obj,
                                thresh_pval = thresh_pval,
                                thresh_ARI = thresh_ARI,
                                thresh_proportion = thresh_proportion)
    valid_multi_k_clust(multi.k.obj)
    return(multi.k.obj)
  } else if (is(x, 'SummExpDR')) {

    if (!is.null(redDimKey)) {
      coords <- t(getEmbeddings(x, key = redDimKey, rows = dims_use, cols = NULL))
      sample.names <- rownames(coords)
      col.names <- colnames(coords)
      if (is.null(dims_use)) {
        dims_use <- 1:ncol(coords)
      }
      # make new data matrix of subsetted PCA embeddings
      coords <- matrix(coords[,dims_use], nrow = nrow(coords), ncol = length(dims_use))
      rownames(coords) <- sample.names
      colnames(coords) <- col.names[dims_use]
      if (is.null(std_norm)) {
        std_norm <- F
      }
    } else {
      coords <- t(SummExpDR::assay(x, i = assay_use))
      if (is.null(std_norm)) {
        std_norm <- T
      }
    }
    # if (is(x@clust.results[[analysis_name]], 'multi_k_clust')) {
    #   warning(paste('overwriting input object\'s pre-existing clustering analysis named', analysis_name))
    # }
    multi.k.obj <- run_multi_k_functions(x = coords,
                                         k_use = k_use,
                                         num_iter = num_iter,
                                         std_norm = std_norm,
                                         algorithm = algorithm,
                                         pct_resample = pct_resample,
                                         affi_K = affi_K,
                                         sigma = sigma,
                                         base_seed = base_seed,
                                         thresh_pval = thresh_pval,
                                         thresh_ARI = thresh_ARI,
                                         thresh_proportion = thresh_proportion,
                                         ...)
    x <- setAnalyses(x, analysis_name, multi.k.obj)
    return(x)
  } else {
    stop(paste('object of class', class(x), 'not supported by this function'))
  }

}

#' Report metrics for quality control, multiple k
#'
#' @inheritParams get_qc_k
#' @details the following QC metrics will be reported:
#' statistical significance for wilcoxon rank sum test comparing pairwise ARI
#' for resampled clustering solutions and those obtained for permuted labels,
#' determined with specified p.value threshold; p.value for said test (is.signif);
#' proportion of pairwise ARI for resampling based clustering solutions with
#' ARI above specified threshold (prop_rs_v_rs);
#' proportion of pairwise ARI in comparisons of reported solution (for given k)
#' to resampling based clusters (prop_rpt_v_rs);
#' does proportion of pairwise ARI values for the reampled clustering comparisons
#' and that calculated for the reported solution x resampled solutions comparisons,
#' exceed a specified threshold proportion (meets_qc);
#' These values are all stored in teh data.frame that is returned
#' @param multi_k_clust multi_k_clust object, output by \code{\link{run_clust_multi_k}}
#' @return multi_k_clust object, with a dataframe with following entries:
#'  - is.signif: logical
#' - prop_rs_v_rs: numeric
#' - prop_rpt_v_rs: numeric
#' - meets_qc: logical
#' - k: integer (# k clusters for partition)
#' see details for more info
#' @export
get_qc_all_k <- function(multi_k_clust,
                         thresh_pval = 0.05,
                         thresh_ARI = 0.90,
                         thresh_proportion = 0.90) {
  values.k <- names(multi_k_clust@single_k_analyses)
  for (i in 1:length(multi_k_clust@single_k_analyses)) {
    k <- values.k[i]
    iter.clust.obj.k <- multi_k_clust@single_k_analyses[[k]]
    qc.rpt.k <- get_qc_k(iter_clust_obj = iter.clust.obj.k,
                         thresh_pval = thresh_pval,
                         thresh_ARI = thresh_ARI,
                         thresh_proportion = thresh_proportion)
    if (i == 1) {
      qc.rpt <- qc.rpt.k
    } else {
      qc.rpt <- rbind(qc.rpt, qc.rpt.k)
    }
  }
  qc.rpt$k <- as.integer(values.k)
  multi_k_clust@report <- qc.rpt
  valid_multi_k_clust(multi_k_clust)
  if (!any(qc.rpt$meets_qc) ) {
    warning('No value k yielded sufficiently stable cluster labels given current thresholds. Consider
            examining the distribution of pairwise ARI among resampled solutions and that between the
            full clustering solution and the resampled solutions')
  }
  return(multi_k_clust)
}

#' Retrieve a clustering solution out of a multi_k_clust object
#'
#' @description retrieve best clustering solution from a given value k,
#' by default, the highest value k that passes 'qc' metrics.
#'
#' @param object multi_k_clust object output by \code{\link{get_qc_all_k}} or
#' \code{\link{run_multi_k_functions}}, or a SummExpDR object output by
#' \code{\link{run_multi_k_functions}}
#' @param pick_best_valid pick best solution for highest value k that passes all qc
#' metrics as described in \code{\link{get_qc_all_k}}
#' @param k NULL, or integer specifying which k to use. if clustering results at
#' given k did not pass qc results, generate a warning but still return result. If
#' specified, return solution reported for clustering analysis with k clusters,
#' @param analysis_use name of clustering analysis, or an integer specifying index
#' of cluster analysis performed. Specifies which multi_k_clust object to retrieve
#' from a SummExpDR object class object
#' @param add_2_metadata add cluster labels to metadata. new metadata column will be
#' named according to the name of the clustering analysis and the value k (# of
#' clusters). E.g. clust.analysis_k3 for clustering analysis 'clust.analysis',
#' best solution for k3
#' (reported by corresponding iter_clust_obj)
#' @return named integer vector of clustering labels, or if object is of class SummExpDR
#' and add_2_metadata = TRUE, a SummExpDR with cluster labels added to metadata
#' @export

get_soln_k <- function(object,
                       pick_best_valid = T,
                       k = NULL,
                       analysis_use = 1,
                       add_2_metadata = F) {

  if (is(object, 'multi_k_clust')) {
    if (pick_best_valid) {
      rpt <- object@report
      if (!(is.data.frame(rpt) && is.logical(rpt$meets_qc))) {
        stop('expected data frame for object@report with columns k and meets_qc.
             was object output by run_multi_k_functions or get_qc_all_k?')
      }
      if (!is.null(k)) {
        warning('ignoring argument k since pick_best_valid set to TRUE')
      }
      if (!any(rpt$meets_qc)) {
        stop('No values k that meet qc thresholds set for this report')
      }
      k_use <- max(rpt$k[which(rpt$meets_qc)])
    } else {
      if (is.numeric(k) && (as.integer(k) == as.numeric(k))) {
        k_use <- k
      } else {
        stop('if pick_best_valid set to false, specify k as an integer value')
      }
    }
    clust.labs <- get_soln_iter_clust(object@single_k_analyses[[as.character(k_use)]],
                                      get_best = T)
    return(clust.labs)
  } else if (is(object, 'SummExpDR')) {
      if (is.numeric(analysis_use) && (analysis_use == as.integer(analysis_use))
          && length(analysis_use) == 1) {
        analyses_keys <- getAnalyses_keys(object)
        analysis_use <- analyses_keys[analysis_use]
      } else if (! (is.character(analysis_use) && (length(analysis_use) == 1)) ) {
        stop('expect analysis_use argument to be either a character of length 1
             or integer lf length 1')
      }
      clust.result <- getAnalyses(object, key = analysis_use)
      clust.labs <- get_soln_k(object = clust.result,
                               pick_best_valid = pick_best_valid,
                               k = k)
      if (add_2_metadata) {
        value.k <- length(unique(clust.labs))
        clust.name <- paste0(analysis_use, '_k', value.k)
        # assume that sample ordering in clusters is same as sample ordering in object
        clust.labs <- as.character(clust.labs)
        names(clust.labs) <- rownames(SummExpDR::colData(object))
        object <- SummExpDR::addColData(object, value = clust.labs, col_name = clust.name)
        return(object)
      } else {
        return(clust.labs)
      }
  }
  # } else if (is(object, 'SummExpDR object')) {
  #   if (is.numeric(analysis_use) && (analysis_use == as.integer(analysis_use))
  #       && length(analysis_use) == 1) {
  #     analysis_use <- names(object@clust.results)[analysis_use]
  #   } else if (! (is.character(analysis_use) && (length(analysis_use) == 1)) ) {
  #     stop('expect analysis_use argument to be either a character of length 1
  #          or integer lf length 1')
  #   }
  #   if (!analysis_use %in% names(object@clust.results)) {
  #     stop(paste('could not find clustering analysis with name', analysis_use))
  #   }
  #   clust.result <- object@clust.results[[analysis_use]]
  #   clust.labs <- get_soln_k(object = clust.result,
  #                            pick_best_valid = pick_best_valid,
  #                            k = k)
  #   if (add_2_metadata) {
  #     value.k <- length(unique(clust.labs))
  #     clust.name <- paste0(analysis_use, '_k', value.k)
  #     col.data <- SummarizedExperiment::colData(object@expt.data)
  #     col.data[,clust.name] <- factor(clust.labs[rownames(col.data)])
  #     SummarizedExperiment::colData(object@expt.data) <- col.data
  #     return(object)
  #   } else {
  #     return(clust.labs)
  #   }
  # } else {
  else {
    stop(paste('object of class', class(object), 'not supported by this function'))
  }
}

