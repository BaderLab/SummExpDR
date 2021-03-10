###############################################################################
### Overlap Utilities ###

#' Fisher's Exact Test, with set input
#'
#' @param set1 character vector containing a set of labels
#' @param set2 character vector containing a set of labels
#' @param master_set set of which set1 and set2 are a subset of and from which they were sampled
#' @description computes fisher's exact test for association between set1 and set2
#' given a master set from which set1 and set 2 were drawn.
#' @return the result of fisher.test run for a contingency table made given the input sets
#' @export
#'
#' @examples
fisher_test_sets <- function(set1, set2, master_set) {
  if (!all(set1 %in% master_set) || !all(set2 %in% master_set)) {
    stop('all members of set1 and set2 must be members of master set')
  }
  TT <- length(intersect(set1, set2))
  TF <- length(set2) - TT
  FT <- length(set1) - TT
  FF <- length(setdiff(master_set, set1)) - TF
  cont.mat <- matrix(c(TT, TF, FT, FF), nrow = 2)
  test.result <- fisher.test(cont.mat, alternative = 'greater')
  return(test.result)
}

#' Calculate overlap between character vectors
#' @param set1 = character vector
#' @param set2 = character vector
#' @param metric = num.ovr for number of overlapping elements in character vectors a nd b,
#' jaccard for jaccard coefficient, ovr.coef for intersect(a,b)/union(a,b),
#' fisher.p for -log10 p.value of fisher's exact test
#' @param master_set = master geneset used for fisher's exact test. set1 and set2 must be subset of this set.
#' only used if fisher.p option selected
#' @return calculated metric of overlap for set1 and set2

set_overlap <- function(set1, set2, metric = 'num.ovr', master_set = NULL) {
  num.common <- length(intersect(set1, set2))
  if (metric == 'num.ovr') {
    return_val <- num.common
  } else if (metric == 'jaccard') {
    num.union <- length(union(set1, set2))
    if (num.union < 1) {
      stop('union of  set1 and set2 is less than 1')
    }
    return_val <- num.common/num.union
  } else if (metric == 'ovr.coef') {
    min.size <- min(length(set1),
                    length(set2))
    if (min.size < 1) {
      stop(paste('expected a minimum set size of 1 for both sets'))
    }
    return_val <- num.common/min.size
  } else if (metric == 'fisher.p') {
    if (!is.character(master_set)) {
      stop('fisher.p specified as method; master_set must be specified as character vector')
    }
    test.result<- fisher_test_sets(set1 = set1,
                                   set2 = set2,
                                   master_set = master_set)
    return_val <- -log10(test.result$p.value)
  } else {
    stop(paste('metric', metric, 'not implemented'))
  }
  return(return_val)
}

#' Calculate Overlap between character all vectors in list(s)
#' @description Given two lists of character vectors, calculate overlap for each
#' combination of character vectors from first and second list, return as matrix. If one list is
#' provided, calculate overlap for each combination of character vectors within given list.
#' @inheritParams set_overlap
#' @param n_cores = number of cores to use in parallel. If NULL, use detected cores on machine - 1.
#' if n_cores set to 1, a serial job is done
#' @param pbar whether or not to use progress bar
#' @return If 2 lists provided, and m x n matrix M where m is number of
#' elements in list1, and n is number of elements in list2, and M[i,j] is
#' calculated overlap between sets list1[[i]] and list2[[j]]. If 1 list provided,
#' an nxn matrix where n is length of list1, M[i,j] is calculated overlap of
#' list1 i and list1 j
#' @export
#'
#' @examples
calc_overlap <- function(list1, list2 = NULL, metric = 'num.ovr', master_set = NULL, n_cores = 1, pbar = TRUE) {

  if (is.null(n_cores)) {
    n_cores <- parallel::detectCores() - 1
  } else {
    tryCatch({
      if (is.numeric(n_cores)) {
        # tolerate floats as long as they denote integer values
        stopifnot(length(n_cores == 1) && as.integer(n_cores) == n_cores)
      } else {
        stop()
      }
    }, error = function(e) {
      stop('n_cores must be integer value or NULL')
    })
  }
  # checks on input
  if (!all(as.logical(lapply(list1, is.character)))) {
    stop('list1 must contain only character vectors')
  } else if (length(list1) < 1) {
    stop('list1 must have nonzero length')
  }

  if (is.null(list2)) {
    list2 <- list1
  } else if (!all(as.logical(lapply(list2, is.character)))) {
    stop('list2 must contain only character vectors')
  } else if (length(list2) < 1) {
    stop('list2 must have nonzero length')
  }

  # loop over lists
  loop_df <- data.frame(tidyr::crossing(data.frame(i = 1:length(list1), stringsAsFactors = FALSE),
                                        data.frame(j = 1:length(list2), stringsAsFactors = FALSE)),
                        stringsAsFactors = FALSE)
  if (n_cores ==  1) {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  } else {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  }
  bp_result <- BiocParallel::bplapply(1:nrow(loop_df), FUN = function(x) {
    i <- loop_df$i[x]
    j <- loop_df$j[x]
    return(set_overlap(list1[[i]], list2[[j]], metric = metric, master_set = master_set))
  }, BPPARAM = bpparam)
  M <- matrix(unlist(bp_result), nrow = length(list1), byrow = TRUE)
  rownames(M) <- names(list1)
  colnames(M) <- names(list2)
  return(M)
}


