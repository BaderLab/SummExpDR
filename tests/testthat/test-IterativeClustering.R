###############################################################################
### Test Clustering Utilities, iterative clustering ###
## note: for whatever reason, need context call at top of file
context('iter_clust_obj definition/validity')
set.seed(5)
DNAm_mat <- t(DNAm_bval)[,1:100]
## helper functions for iter_clust_obj methods
context('make_inp_mat')
test_that('make_inp_mat works properly, snf_affi mode', {
  for (std_norm in c(T,F)) {
    if (std_norm) {
      ref_norm <- SNFtool::standardNormalization(DNAm_mat)
    } else {
      ref_norm <- DNAm_mat
    }
    for (as_kernMat in c(T,F)) {
      ref_dist <- SNFtool::dist2(ref_norm, ref_norm)^(1/2)
      ref_mat <- SNFtool::affinityMatrix(ref_dist, K = 20, sigma = 0.5)
      test_mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = 'snf_affi', affi_K = 20, sigma = 0.5,
                                 as_kernMat = as_kernMat)
      if (as_kernMat) {
        ref_mat <- kernlab::as.kernelMatrix(ref_mat)
      }
      expect_identical(test_mat, ref_mat)
      expect_identical(rownames(test_mat), rownames(DNAm_mat))
      expect_identical(colnames(test_mat), rownames(DNAm_mat))
    }
  }
})
test_that('make_inp_mat works properly, mxn mode', {

  for (std_norm in c(T,F)) {
    if (std_norm) {
      ref_norm <- SNFtool::standardNormalization(DNAm_mat)
    } else {
      ref_norm <- DNAm_mat
    }

    test_mat <- make_input_mat(DNAm_mat,
                               std_norm = std_norm,
                               affi_K = 20,
                               sigma = 0.5,
                               mode = 'mxn')
    expect_identical(test_mat, ref_norm)
    expect_identical(rownames(test_mat), rownames(DNAm_mat))
  }
})

context('clust_once')
test_that('clust_once works properly', {

  seed <- 2

  for (std_norm in c(T,F)) {
    for (k in 2:5) {

      mxn.mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = 'mxn')
      snf.affi.mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = 'snf_affi', affi_K = 20, sigma = 0.5)
      snf.affi.kernlab.mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = 'snf_affi',
                                             affi_K = 20, sigma = 0.5, as_kernMat = T)

      set.seed(seed)
      kmeans.ref <- kmeans(x = mxn.mat, centers = k)$cluster
      kmeans.res <- clust_once(inp_mat = mxn.mat,
                               algorithm = 'kmeans',
                               k = k,
                               seed)
      expect_identical(kmeans.res, kmeans.ref)

      spectral.snf.ref <- SNFtool::spectralClustering(snf.affi.mat, K = k, type = 3)
      spectral.snf.res <- clust_once(snf.affi.mat, algorithm = 'spectral_SNFtool', k = k, seed)
      expect_identical(spectral.snf.res, spectral.snf.ref)

      set.seed(seed)
      spectral.kernlab.ref <- kernlab::specc(x = snf.affi.kernlab.mat, centers = k)@.Data
      spectral.kernlab.res <- clust_once(inp_mat = snf.affi.kernlab.mat,
                                         algorithm = 'spectral_kernlab',
                                         k = k,
                                         seed)
      expect_identical(spectral.kernlab.res, spectral.kernlab.ref)

    }
  }

})

context('iter_clust function')
test_that('iter_clust produces correct output, spectral_SNFtool mode, no resampling', {
  ## can easily change the name of this test to 'non-restart' modes
  ## such as hierarchical clustering
  for (std_norm in c(T,F)) {
    k <- 10
    affi_K <- 20
    sigma <- 0.5
    ## note that iterations here is specified but should not
    ## result in a 10xnum_sample matrix as spectral mode should
    ## only perform clustering once
    iterations <- 10
    inp_mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = 'snf_affi',
                              affi_K = affi_K, sigma = sigma)
    ref_clust <- SNFtool::spectralClustering(affinity = inp_mat, K = k, type = 3)
    names(ref_clust) <- rownames(DNAm_mat)
    test_clust <- iter_clust(DNAm_mat, num_iter = iterations,
                             algorithm = 'spectral_SNFtool',
                             std_norm = std_norm,
                             k = k,
                             affi_K = affi_K,
                             sigma = sigma,
                             resample = F)
    expect_equal(colnames(test_clust), names(ref_clust))
    expect_equal(test_clust[1,], ref_clust)
    expect_equal(nrow(test_clust), 1)
    expect_equal(class(test_clust), 'matrix')
  }

})
test_that('iter_clust produces correct output, restart modes', {

  algorithms.use <- c('spectral_kernlab',  'kmeans')

  ## restart modes include kmeans and spectral_kernlab
  iterations <- 10
  k <- 3
  base_seed <- 4
  sample_names <- rownames(DNAm_mat)
  num_samples <- length(sample_names)
  for (algorithm in algorithms.use) {
    for (std_norm in c(T,F)) {

      ## switch input matrix type. can change to 'algorithm %in% c(algorithm1, algorithm2)
      ## for algorithms that require identical input type
      if (algorithm == 'kmeans')  {
        inp.mat.type <- 'mxn'
        as_kernMat <- F
      } else if (algorithm == 'spectral_kernlab') {
        inp.mat.type <- 'snf_affi'
        as_kernMat <- T
      }
      inp_mat <- make_input_mat(DNAm_mat, std_norm = std_norm, mode = inp.mat.type, as_kernMat = as_kernMat)
      test_mat <- iter_clust(DNAm_mat,
                             num_iter = iterations,
                             std_norm = std_norm,
                             algorithm = algorithm,
                             k = k,
                             base_seed = base_seed,
                             resample = F)
      ## test dimensions/colnames
      expect_equal(colnames(test_mat), sample_names)
      expect_equal(dim(test_mat), c(iterations, num_samples))
      ## test comments
      for (i in 1:nrow(test_mat)) {
        seed.use <- base_seed + i - 1
        clust.res <- clust_once(inp_mat = inp_mat,
                                algorithm = algorithm,
                                k = k,
                                seed = seed.use)
        expect_true(all(test_mat[i,] == clust.res))
      }
    }
  }
})
test_that('iter_clust produces correct output, any mode, resampling', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (algorithm in algorithms.use) {
    for (std_norm in c(T,F)) {
      iterations <- 10
      k <- 10
      affi_K <- 20
      sigma = 0.5
      base_seed <- 4
      sample_names <- rownames(DNAm_mat)
      num_samples <- length(sample_names)
      pct.resample <- 0.70
      samples_keep <- floor(num_samples*pct.resample)
      test_mat <- iter_clust(DNAm_mat,
                             num_iter = iterations,
                             algorithm = algorithm,
                             base_seed = base_seed,
                             k = k,
                             affi_K = affi_K,
                             sigma = sigma,
                             resample = T,
                             pct_resample = pct.resample,
                             std_norm = std_norm)
      ## test dimensions/colnames
      expect_equal(colnames(test_mat), sample_names)
      expect_equal(dim(test_mat), c(iterations, num_samples))
      ## test contents
      for (i in 1:iterations) {
        ## expect number of unique elements in test_mat to be equal to
        ## k + 1 (k clusters + 0 label). Note that due to implementation
        ## of spectral clustering from SNFtool, possible to get less than
        ## k clusters
        # expect_equal(length(unique(test_mat[i,])), k + 1)
        expect_true(all(unique(test_mat[i,]) %in% c(1:k, 0)))
        ## check number of nonzero entries against number of samples to keep
        expect_equal(sum(as.numeric(!test_mat[i,] == 0)), samples_keep)
        ## expect different result in each iteration
        if (i > 1) {
          expect_false(all(test_mat[i,] == test_mat[i - 1,]))
        }
        ## check one row's entries. Are cluster labels for given subset identical in reference + function output?
        seed.use <- base_seed + i - 1
        set.seed(seed.use)
        samples_subs <- sample_names[sample(1:num_samples, samples_keep, replace = F)]
        as_kernMat <- ifelse((algorithm == 'spectral_kernlab'), T, F)
        inp.mat.type <- ifelse((algorithm %in% c('spectral_kernlab', 'spectral_SNFtool')), 'snf_affi', 'mxn')
        inp_mat <- make_input_mat(DNAm_mat[samples_subs,], std_norm = std_norm,
                                  mode = inp.mat.type, as_kernMat = as_kernMat)
        clust.ref <- clust_once(inp_mat = inp_mat,
                                algorithm = algorithm,
                                k = k,
                                seed = seed.use)
        expect_true(all(clust.ref == test_mat[i, samples_subs]))
        ## are all other entries 0?
        expect_true(all(test_mat[i, !sample_names %in% samples_subs] == 0))
      }
    }
  }
})

test_that('iter_clust produces different labels with different k, all modes', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (clust_mode in algorithms.use) {
    for (resample in c(T, F)) {
      for (std_norm in c(T,F))
        k10_result <- iter_clust(DNAm_mat,
                                 k = 10,
                                 base_seed = 12345,
                                 std_norm = std_norm,
                                 resample = resample,
                                 algorithm = clust_mode,
                                 num_iter = 10)
      k11_result <- iter_clust(DNAm_mat,
                               k = 11,
                               std_norm = std_norm,
                               base_seed = 12345,
                               resample = resample,
                               algorithm = clust_mode,
                               num_iter = 10)
      expect_false(identical(k10_result, k11_result))
      ## check number of unique elments in each matrix
      expect_false(length(unique(as.numeric(k10_result)))
                   == length(unique(as.numeric(k11_result))))
    }
  }
})

context('perm_labels')
test_that('perm_labels produces expected output', {
  for (resample in c(T)) {
    base_seed = 2
    clust_result <- iter_clust(DNAm_mat,
                               k = 3,
                               num_iter = 10,
                               resample = resample,
                               base_seed = base_seed)
    ref_perm <- clust_result
    for (i in 1:nrow(clust_result)) {
      set.seed((base_seed + i - 1)^2)
      ref_perm[i,] <- gtools::permute(ref_perm[i,])
    }
    test_perm <- perm_labels(clust_result, base_seed = 2)
    expect_identical(test_perm, ref_perm)
    expect_equal(length(unique(test_perm)), length(unique(clust_result)))
  }
})

test_that('perm_labels output determined by seed', {
  for (resample in c(T)) {
    base_seed = 2
    clust_result <- iter_clust(DNAm_mat,
                               k = 3,
                               num_iter = 10,
                               resample = resample,
                               base_seed = base_seed)
    perm_a  <- perm_labels(M = clust_result, base_seed = 12345)
    perm_b <- perm_labels(M = clust_result, base_seed = 2)
    perm_c <- perm_labels(M = clust_result, base_seed = 12345)
    expect_identical(perm_a, perm_c)
    expect_false(identical(perm_a, perm_b))
  }
})
context('get_sim_mat')
test_that('get_sim_mat works properly', {
  ## note that we only test with kmeans mode as all we need
  ## is cluster labels to be spit out
  k = 10
  iterations = 10
  base_seed = 2
  pct.resample = 0.80
  ## test case of no resampling
  clust.res.full <- iter_clust(x = DNAm_mat,
                               k = k,
                               algorithm = 'kmeans',
                               num_iter = iterations,
                               resample = F,
                               pct_resample = pct.resample,
                               base_seed = base_seed)
  sim.full.data <- get_sim_mat(clust.res.full)
  ARI.3.4 <- mclust::adjustedRandIndex(clust.res.full[3,], clust.res.full[4,])
  expect_equal(sim.full.data[3,4], ARI.3.4)
  ## test case of resampling
  clust.res.rs <- iter_clust(x = DNAm_mat,
                             k = k,
                             algorithm = 'kmeans',
                             num_iter = iterations,
                             resample = T,
                             pct_resample = pct.resample,
                             base_seed = base_seed)
  sim.rs.data <- get_sim_mat(clust.res.rs)
  common.5.7 <- which((clust.res.rs[5,] != 0) & (clust.res.rs[7, ] != 0))
  ARI.5.7 <- mclust::adjustedRandIndex(clust.res.rs[5, common.5.7], clust.res.rs[7, common.5.7])
  expect_equal(ARI.5.7, sim.rs.data[5,7])
  ## check diagonal of each ARI matrix
  for (i in 1:iterations) {
    expect_equal(sim.full.data[i,i], 1)
    expect_equal(sim.rs.data[i,i], 1)
  }

})

## object methods
context('iter_clust_obj methods')
test_that('do_clust works properly', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (clust_mode in algorithms.use) {
    for (std_norm in c(T,F)) {
      k = 10
      affi_K = 20
      sigma = 0.5
      iterations = 10
      base_seed = 2
      pct.resample = 0.80

      new_obj <- make_iter_clust_obj()
      new_obj <- do_clust(object = new_obj,
                          inp_data = DNAm_mat,
                          std_norm = std_norm,
                          k = k,
                          algorithm = clust_mode,
                          affi_K = affi_K,
                          sigma = sigma,
                          num_iter = iterations,
                          pct_resample = pct.resample,
                          base_seed = base_seed)
      ref_clust_full <- iter_clust(x = DNAm_mat,
                                   k = k,
                                   algorithm = clust_mode,
                                   affi_K = affi_K,
                                   sigma = sigma,
                                   num_iter = iterations,
                                   std_norm = std_norm,
                                   resample = F,
                                   pct_resample = pct.resample,
                                   base_seed = base_seed)
      ref_clust_rs <- iter_clust(x = DNAm_mat,
                                 k = k,
                                 algorithm = clust_mode,
                                 affi_K = affi_K,
                                 sigma = sigma,
                                 num_iter = iterations,
                                 std_norm = std_norm,
                                 resample = T,
                                 pct_resample = pct.resample,
                                 base_seed = base_seed)
      ref_clust_rs_perm <- perm_labels(ref_clust_rs,
                                       base_seed = base_seed)
      expect_identical(new_obj@clusts_full, ref_clust_full)
      expect_identical(new_obj@clusts_rs, ref_clust_rs)
      expect_identical(new_obj@clusts_rs_perm, ref_clust_rs_perm)
    }
  }
})

test_that('calc_sim_mats works properly', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (clust_mode in algorithms.use) {
    k = 10
    iterations = 10
    base_seed = 2
    pct.resample = 0.80

    new_obj <- make_iter_clust_obj()
    new_obj <- do_clust(object = new_obj,
                        inp_data = DNAm_mat,
                        k = k,
                        algorithm = clust_mode,
                        num_iter = iterations,
                        pct_resample = pct.resample,
                        base_seed = base_seed)
    new_obj_ref <- new_obj
    new_obj <- calc_sim_mats(new_obj)

    ## if a restart clustering mode, sim_mat_full should have
    ## a matrix value equivalent to output of get_sim_mat.
    ## otherwise should be NULL
    if (clust_mode %in% c('spectral_kernlab',  'kmeans')) {
      new_obj_ref@sim_mat_full <- get_sim_mat(new_obj_ref@clusts_full)
      expect_identical(new_obj@sim_mat_full, new_obj_ref@sim_mat_full)
    } else {
      expect_true(is.null(new_obj@sim_mat_full))
    }

    new_obj_ref@sim_mat_rs <- get_sim_mat(new_obj_ref@clusts_rs)
    new_obj_ref@sim_mat_rs_perm <- get_sim_mat(new_obj_ref@clusts_rs_perm)
    expect_identical(new_obj@sim_mat_rs, new_obj_ref@sim_mat_rs)
    expect_identical(new_obj@sim_mat_rs_perm, new_obj_ref@sim_mat_rs_perm)
  }
})

test_that('unique_solutions works properly', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (clust_mode in algorithms.use) {
    k <- 10
    iterations <- 50
    base_seed <- 2
    pct.resample <- 0.80
    new_obj <- make_iter_clust_obj()
    new_obj <- do_clust(object = new_obj,
                        inp_data = DNAm_mat,
                        k = k,
                        algorithm = clust_mode,
                        num_iter = iterations,
                        pct_resample = pct.resample,
                        base_seed = base_seed)
    new_obj <- calc_sim_mats(new_obj)
    new_obj <- unique_solutions(new_obj)
    unique_solns <- new_obj@unique_solns
    soln_metadata <- new_obj@soln_metadata
    expect_true(is.matrix(unique_solns))
    expect_true(is.data.frame(soln_metadata))
    if (clust_mode == 'spectral_SNFtool') {
      expect_equal(nrow(unique_solns), 1)
      expect_equal(nrow(soln_metadata), 1)
      expect_equal(soln_metadata[1,'solution'], 1)
      expect_equal(soln_metadata[1,'proportion'], 1)
    } else {
      ## expect each solution found in iterations to be identical
      ## to 1 and only 1 of the unique solutions. Also count number of
      ## members per each solution
      num.solns <- numeric(nrow(unique_solns))
      for (i in 1:iterations) {
        clust.soln.i <- new_obj@clusts_full[i,]
        num.matches <- 0
        for (j in 1:nrow(unique_solns)) {
          unique_solns.j <- unique_solns[j,]
          if (mclust::adjustedRandIndex(clust.soln.i, unique_solns.j) == 1) {
            num.matches <- num.matches + 1
            num.solns[j] <- num.solns[j] + 1
          }

        }
        if (!num.matches == 1) {
          stop(paste0('Expected number of unique solutnion matches for a given
                      clustering solution to be 1 and only 1, at clustering ',
                      i, ' got ', num.matches, 'matches'))
        }
      }
      ## calculated proportion of all clusterings, for each unique solution
      ## found
      ref.prop <- num.solns/iterations
      ## run expectations
      expect_identical(soln_metadata$solution, 1:nrow(unique_solns))
      expect_equal(soln_metadata$proportion, ref.prop)
      expect_true(all(soln_metadata$proportion > 0))
      expect_true(all(soln_metadata$proportion < 1))
      expect_true(sum(soln_metadata$proportion) == 1)
    }
  }
})

test_that('compare_rs_perm passes sanity check', {
  ## run clustering with well separated data (2 gaussian distributions).
  ## should have very similar resample based clustering results
  iterations <- 50
  base_seed <- 2
  pct.resample <- 0.80

  gauss_2 <- DNAm_mat
  num_samples <- nrow(DNAm_mat)
  first_half <- 1:(floor(num_samples/2))
  set.seed(base_seed)
  gauss_2[first_half,] <- apply(gauss_2[first_half,],
                                MARGIN = 2,
                                FUN = function(x) {
                                  return(rnorm(length(x),
                                               mean = 3,
                                               sd = 1))
                                })
  second_half <- ceiling(num_samples/2):num_samples
  set.seed(base_seed)
  gauss_2[second_half,] <- apply(gauss_2[second_half,],
                                 MARGIN = 2,
                                 FUN = function(x) {
                                   return(rnorm(length(x),
                                                mean = -3,
                                                sd = 1))
                                 })

  algorithms.use <- c('kmeans')

  for (algorithm in algorithms.use) {
    gauss_obj <- make_iter_clust_obj()
    gauss_obj <- do_clust(object = gauss_obj,
                          inp_data = gauss_2,
                          k = 2,
                          algorithm = algorithm,
                          num_iter = iterations,
                          pct_resample = pct.resample,
                          base_seed = base_seed)
    gauss_obj <- calc_sim_mats(gauss_obj)
    gauss_obj <- compare_rs_perm(gauss_obj)
    wcox_result <- gauss_obj@sim_rs_v_perm
    expect_is(wcox_result, 'htest')
    expect_lt(wcox_result$p.value, 0.05)
  }

})

## For now skipping tests of some of simpler functions + methods

test_that('run_iter_clust_functions runs without error,
          produces expected output, reports correct solutions', {

            algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

            for (clust_mode in algorithms.use) {
              output_obj <- run_iter_clust_functions(inp_data = DNAm_mat,
                                                     k = 10,
                                                     num_iter = 50,
                                                     algorithm = clust_mode,
                                                     pct_resample = 0.80,
                                                     affi_K = 20,
                                                     sigma = 0.5,
                                                     base_seed = 2)
              ## sanity checks on reported solution + solution metadata
              soln_metadata <- get_soln_metadata(output_obj)
              expect_equal(soln_metadata$med_ARI_rs, apply(output_obj@compar_mat_full_rs,
                                                           MARGIN = 1,
                                                           median))
              soln.ind <- which(soln_metadata$soln_report)
              unique_solns <- output_obj@unique_solns
              expect_equal(length(soln.ind), 1)
              expect_equal(soln_metadata[soln.ind, 'med_ARI_rs'], max(soln_metadata$med_ARI_rs))
              max.med.ARI.solns <- which(soln_metadata$med_ARI_rs ==  max(soln_metadata$med_ARI_rs))
              expect_equal(soln_metadata$proportion[soln.ind],
                           max(soln_metadata[max.med.ARI.solns, 'proportion']))
              soln.retrieved <- get_soln_iter_clust(output_obj, get_best = T)
              expect_equal(soln.retrieved, unique_solns[soln.ind,])
              expect_equal(names(soln.retrieved), colnames(unique_solns))
            }
          })

context('quality control')
test_that('get_qc_k sanity check', {

  algorithms.use <- c('spectral_SNFtool', 'spectral_kernlab',  'kmeans')

  for (clust_mode in algorithms.use) {
    non.sig.cts <- 0
    sig.cts <- 0
    for (k in 2:5) {
      output.obj <- run_iter_clust_functions(inp_data = DNAm_mat,
                                             k = k,
                                             num_iter = 50,
                                             algorithm = clust_mode,
                                             pct_resample = 0.80,
                                             affi_K = 20,
                                             sigma = 0.5,
                                             base_seed = 2)
      result <- get_qc_k(output.obj,
                         thresh_pval = 0.05,
                         thresh_ARI = 0.6,
                         thresh_proportion = 0.5)
      sim_dist_rs_v_rs <- get_sim_dist(output.obj@sim_mat_rs)
      sim_dist_rpt_v_rs <- as.numeric(output.obj@compar_mat_full_rs)
      expect_equal(result$prop_rs_v_rs,
                   sum(sim_dist_rs_v_rs > 0.6)/length(sim_dist_rs_v_rs))
      expect_equal(result$prop_rpt_v_rs,
                   sum(sim_dist_rpt_v_rs > 0.6)/length(sim_dist_rpt_v_rs))
      if (result$meets_qc) {
        expect_true(result$is_signif)
        expect_lt(output.obj@sim_rs_v_perm$p.value,
                  0.05)
        expect_gt(result$prop_rs_v_rs, 0.5)
        expect_gt(result$prop_rpt_v_rs, 0.5)
        sig.cts <- sig.cts + 1
      } else {
        non.sig.cts <- non.sig.cts + 1
      }
    }
    if (sig.cts <= 0) {
      warning('expected at least 1 test on a passing case for get_qc_k')
    }
    if (non.sig.cts <= 0) {
      warning('expected at least 1 test on a non-passing case for get_qc_k')
    }
  }

})


context('multi_k_clust_obj class definition')

test_that('multi_k_clust_obj vallidity function',  {
  multi.k <- create_multi_k_clust()
  bad.single_k_analyses <- list('asdf', 2.5, 1:5, list('asdf','ghjk'),
                                list(1,2,3), list(cars, mtcars))
  for (element in bad.single_k_analyses) {
    multi.k@single_k_analyses <- element
    expect_error(valid_multi_k_clust(multi.k))
  }
  bad.reports <- bad.single_k_analyses
  for (element in bad.reports) {
    multi.k@report <- element
    expect_error(valid_multi_k_clust(multi.k))
  }
})

context('multi_k_clust_obj clustering')

## Setup data for this context and context "get_soln_k works properly"
DNAm_SE <- SummarizedExperiment::SummarizedExperiment(assays = list(bval = t(DNAm_mat)))
test_obj <- create_SummExpDR(DNAm_SE)
test_obj <- runPCA(test_obj, std_norm = TRUE, i = 'bval')
pca_coords <- as.matrix(fetchData(test_obj, redDimKeys = c('PCA'), varsFetch = paste0('PC', 1:10)))

## clustering algorithms to use
algorithms.use <- c('spectral_kernlab', 'spectral_SNFtool', 'kmeans')
k_use <- 2:3
num_iter <- 25
pct_resample <- 0.80
affi_K <- 8
sigma <- 0.5
thresh_pval <- 0.10
thresh_ARI <- 0.50
thresh_proportion <- 0.50
base_seed <- 2

test_that('run_clust_multi_k produces expected clustering output', {
  for (algorithm in algorithms.use) {
    for (std_norm in c(T,F)) {
      multi.k <- create_multi_k_clust()
      multi.k <- run_clust_multi_k(object = multi.k,
                                   inp_data = DNAm_mat,
                                   k_use = k_use,
                                   std_norm = std_norm,
                                   num_iter = num_iter,
                                   algorithm = algorithm,
                                   pct_resample = pct_resample,
                                   affi_K = affi_K,
                                   sigma = sigma,
                                   base_seed = base_seed)
      for (i in k_use) {
        ref.i <- run_iter_clust_functions(inp_data = DNAm_mat,
                                          k = i,
                                          std_norm = std_norm,
                                          num_iter = num_iter,
                                          algorithm = algorithm,
                                          pct_resample = pct_resample,
                                          affi_K = affi_K,
                                          sigma = sigma,
                                          base_seed = base_seed)
        expect_equal(ref.i, multi.k@single_k_analyses[[as.character(i)]])
      }
    }
  }
})

test_that('run_multi_k_functions produces expected clustering output', {

  for (algorithm in algorithms.use) {
    for (std_norm in c(T,F)) {
      ## test clustering on data matrix
      multi.k.ref <- create_multi_k_clust()
      multi.k.ref <- run_clust_multi_k(object = multi.k.ref,
                                       inp_data = DNAm_mat,
                                       std_norm = std_norm,
                                       k_use = k_use,
                                       num_iter = num_iter,
                                       algorithm = algorithm,
                                       pct_resample = pct_resample,
                                       affi_K = affi_K,
                                       sigma = sigma,
                                       base_seed = base_seed)
      multi.k.ref <- get_qc_all_k(multi_k_clust = multi.k.ref,
                                  thresh_pval = thresh_pval,
                                  thresh_ARI = thresh_ARI,
                                  thresh_proportion = thresh_proportion)
      multi.k.mat <- run_multi_k_functions(x = DNAm_mat,
                                           std_norm = std_norm,
                                           analysis_name = 'clust.analysis',
                                           k_use = k_use,
                                           num_iter = num_iter,
                                           algorithm = algorithm,
                                           pct_resample = pct_resample,
                                           affi_K = affi_K,
                                           sigma = sigma,
                                           base_seed = base_seed,
                                           thresh_pval = thresh_pval,
                                           thresh_ARI = thresh_ARI,
                                           thresh_proportion = thresh_proportion)

      expect_identical(multi.k.mat, multi.k.ref)

      ## test clustering using SummExpDR input
      multi.k.SummExpDR.data <- run_multi_k_functions(x = test_obj,
                                                      analysis_name = 'clust.analysis',
                                                      redDimKey = NULL,
                                                      dims_use = NULL,
                                                      std_norm = std_norm,
                                                      k_use = k_use,
                                                      num_iter = num_iter,
                                                      algorithm = algorithm,
                                                      pct_resample = pct_resample,
                                                      affi_K = affi_K,
                                                      sigma = sigma,
                                                      base_seed = base_seed,
                                                      thresh_pval = thresh_pval,
                                                      thresh_ARI = thresh_ARI,
                                                      thresh_proportion = thresh_proportion)
      clust_result <- getAnalyses(multi.k.SummExpDR.data, 'clust.analysis')
      expect_identical(clust_result, multi.k.ref)

      ## test clustering on principal components
      multi.k.ref.pca <- create_multi_k_clust()
      multi.k.ref.pca <- run_clust_multi_k(object = multi.k.ref.pca,
                                           inp_data = pca_coords[,1:10],
                                           std_norm = std_norm,
                                           k_use = k_use,
                                           num_iter = num_iter,
                                           algorithm = algorithm,
                                           pct_resample = pct_resample,
                                           affi_K = affi_K,
                                           sigma = sigma,
                                           base_seed = base_seed)
      multi.k.ref.pca <- get_qc_all_k(multi.k.ref.pca,
                                      thresh_pval = thresh_pval,
                                      thresh_ARI = thresh_ARI,
                                      thresh_proportion = thresh_proportion)
      multi.k.SummExpDR.pca <- run_multi_k_functions(x = test_obj,
                                                     analysis_name = 'clust.analysis',
                                                     redDimKey = 'PCA',
                                                     dims_use = 1:10,
                                                     std_norm = std_norm,
                                                     k_use = k_use,
                                                     num_iter = num_iter,
                                                     algorithm = algorithm,
                                                     pct_resample = pct_resample,
                                                     affi_K = affi_K,
                                                     sigma = sigma,
                                                     base_seed = base_seed,
                                                     thresh_pval = thresh_pval,
                                                     thresh_ARI = thresh_ARI,
                                                     thresh_proportion = thresh_proportion)
      clust_result <- getAnalyses(multi.k.SummExpDR.pca, 'clust.analysis')
      expect_identical(clust_result, multi.k.ref.pca)
    }
  }

})

# test_that('default behavior of run_multi_k_functions matches expectation', {
#
#   for (algorithm in algorithms.use) {
#
#     ## test clustering on data matrix
#     multi.k.ref <- create_multi_k_clust()
#     multi.k.ref <- run_clust_multi_k(object = multi.k.ref,
#                                      inp_data = DNAm_mat,
#                                      std_norm = T,
#                                      k_use = k_use,
#                                      num_iter = num_iter,
#                                      algorithm = algorithm,
#                                      pct_resample = pct_resample,
#                                      affi_K = affi_K,
#                                      sigma = sigma,
#                                      base_seed = base_seed)
#     multi.k.ref <- get_qc_all_k(multi_k_clust = multi.k.ref,
#                                 thresh_pval = thresh_pval,
#                                 thresh_ARI = thresh_ARI,
#                                 thresh_proportion = thresh_proportion)
#     multi.k.mat <- run_multi_k_functions(x = DNAm_mat,
#                                          std_norm = NULL,
#                                          analysis_name = 'clust.analysis',
#                                          k_use = k_use,
#                                          num_iter = num_iter,
#                                          algorithm = algorithm,
#                                          pct_resample = pct_resample,
#                                          affi_K = affi_K,
#                                          sigma = sigma,
#                                          base_seed = base_seed,
#                                          thresh_pval = thresh_pval,
#                                          thresh_ARI = thresh_ARI,
#                                          thresh_proportion = thresh_proportion)
#
#     expect_identical(multi.k.mat, multi.k.ref)
#
#     ## test clustering using SummExpDR input
#     multi.k.SummExpDR.data <- run_multi_k_functions(x = test_obj,
#                                                   analysis_name = 'clust.analysis',
#                                                   redDimKey = NULL,
#                                                   dims_use = NULL,
#                                                   std_norm = NULL,
#                                                   k_use = k_use,
#                                                   num_iter = num_iter,
#                                                   algorithm = algorithm,
#                                                   pct_resample = pct_resample,
#                                                   affi_K = affi_K,
#                                                   sigma = sigma,
#                                                   base_seed = base_seed,
#                                                   thresh_pval = thresh_pval,
#                                                   thresh_ARI = thresh_ARI,
#                                                   thresh_proportion = thresh_proportion)
#     clust_result <- getAnalyses(multi.k.SummExpDR.data, 'clust.analysis')
#     expect_identical(clust_result, multi.k.ref)
#
#     ## test clustering on principal components
#     multi.k.ref.pca <- create_multi_k_clust()
#     multi.k.ref.pca <- run_clust_multi_k(object = multi.k.ref.pca,
#                                          inp_data = pca_coords[,1:10],
#                                          std_norm = F,
#                                          k_use = k_use,
#                                          num_iter = num_iter,
#                                          algorithm = algorithm,
#                                          pct_resample = pct_resample,
#                                          affi_K = affi_K,
#                                          sigma = sigma,
#                                          base_seed = base_seed)
#     multi.k.ref.pca <- get_qc_all_k(multi.k.ref.pca,
#                                     thresh_pval = thresh_pval,
#                                     thresh_ARI = thresh_ARI,
#                                     thresh_proportion = thresh_proportion)
#     multi.k.SummExpDR.pca <- run_multi_k_functions(x = test_obj,
#                                                  analysis_name = 'clust.analysis',
#                                                  redDimKey = 'PCA',
#                                                  dims_use = 1:10,
#                                                  std_norm = NULL,
#                                                  k_use = k_use,
#                                                  num_iter = num_iter,
#                                                  algorithm = algorithm,
#                                                  pct_resample = pct_resample,
#                                                  affi_K = affi_K,
#                                                  sigma = sigma,
#                                                  base_seed = base_seed,
#                                                  thresh_pval = thresh_pval,
#                                                  thresh_ARI = thresh_ARI,
#                                                  thresh_proportion = thresh_proportion)
#     clust_result <- getAnalyses(multi.k.SummExpDR.pca, 'clust.analysis')
#     expect_identical(clust_result, multi.k.ref.pca)
#
#   }
# })

context('get_soln_k')
test_that('get_soln_k functions properly', {

  algorithm <- 'kmeans'

  for (k in 1:2) {

    ## try with both specified k and with NULL value
    if (k == 1) {
      k <- NULL
    }

    multi.k.clust <- run_multi_k_functions(x = DNAm_mat,
                                           analysis_name = 'clust.analysis',
                                           k_use = k_use,
                                           num_iter = num_iter,
                                           algorithm = algorithm,
                                           pct_resample = pct_resample,
                                           affi_K = affi_K,
                                           sigma = sigma,
                                           base_seed = base_seed,
                                           thresh_pval = thresh_pval,
                                           thresh_ARI = thresh_ARI,
                                           thresh_proportion = thresh_proportion)

    if (is.null(k)) {
      pick_best_valid <- T
    } else {
      pick_best_valid <- F
    }
    soln.multi.k.clust <- get_soln_k(object = multi.k.clust,
                                     pick_best_valid = pick_best_valid,
                                     k = k)
    ## testing that solution output for multi.k.clust input is correct
    k.found <- length(unique(soln.multi.k.clust))
    if (pick_best_valid) {
      expect_equal(k.found,
                   multi.k.clust@report[ max(which(multi.k.clust@report$meets_qc)), 'k'])
    } else {
      expect_equal(k.found, k)
    }

    expect_equal(soln.multi.k.clust,
                 get_soln_iter_clust(multi.k.clust@single_k_analyses[[as.character(k.found)]],
                                     get_best = T))

    SummExpDR.obj <- run_multi_k_functions(x = test_obj,
                                          analysis_name = 'clust.analysis',
                                          redDimKey = NULL,
                                          dims_use = NULL,
                                          k_use = k_use,
                                          num_iter = num_iter,
                                          algorithm = algorithm,
                                          pct_resample = pct_resample,
                                          affi_K = affi_K,
                                          sigma = sigma,
                                          base_seed = base_seed,
                                          thresh_pval = thresh_pval,
                                          thresh_ARI = thresh_ARI,
                                          thresh_proportion = thresh_proportion)

    soln.SummExpDR.obj <- get_soln_k(object = SummExpDR.obj,
                                pick_best_valid = pick_best_valid,
                                k = k,
                                add_2_metadata = F)


    ## test that cluster labels match output expected of multi_k_clust input
    expect_equal(soln.SummExpDR.obj, soln.multi.k.clust)

    SummExpDR.obj.w.meta <- get_soln_k(object = SummExpDR.obj,
                                  pick_best_valid = pick_best_valid,
                                  k = k,
                                  add_2_metadata = T)
    ## test that cluster labels properly assigned to metadata if add_2_metadata = T
    col.data <- SummarizedExperiment::colData(getSummExp(SummExpDR.obj.w.meta))
    expect_equal(as.vector(col.data[, grep('clust.analysis', colnames(col.data))]),
                 as.character(soln.SummExpDR.obj[rownames(col.data)]))

    ## results of specifying analysis_use correctly should be identical to the
    ## output using first clustering analysis in clust.results slot
    expect_equal(SummExpDR.obj.w.meta, get_soln_k(object = SummExpDR.obj,
                                                 pick_best_valid = pick_best_valid,
                                                 k = k,
                                                 add_2_metadata = T,
                                                 analysis_use = 'clust.analysis'))
    expect_equal(soln.SummExpDR.obj, get_soln_k(object = SummExpDR.obj,
                                               pick_best_valid = pick_best_valid,
                                               k = k,
                                               add_2_metadata = F,
                                               analysis_use = 'clust.analysis'))

  }

})

test_that('get_soln_k fails on unexpected input or expected failure condition', {
  algorithm <- 'kmeans'

  SummExpDR.obj <- run_multi_k_functions(x = test_obj,
                                        analysis_name = 'clust.analysis',
                                        redDimKey = NULL,
                                        dims_use = NULL,
                                        k_use = k_use,
                                        num_iter = num_iter,
                                        algorithm = algorithm,
                                        pct_resample = pct_resample,
                                        affi_K = affi_K,
                                        sigma = sigma,
                                        base_seed = base_seed,
                                        thresh_pval = thresh_pval,
                                        thresh_ARI = thresh_ARI,
                                        thresh_proportion = thresh_proportion)

  ## test expected failures based on non-existent clustering analyses
  expect_error(get_soln_k(object = SummExpDR.obj,
                          pick_best_valid = F,
                          k = 2,
                          add_2_metadata = F,
                          analysis_use = 2))
  expect_error(get_soln_k(object = SummExpDR.obj,
                          pick_best_valid = F,
                          k = 2,
                          add_2_metadata = F,
                          analysis_use = 'asdf'))

  ## test for failure on object with no qc report

  expect_error(get_soln_k(test_obj))

  ## test for failure on multi_k_clust_obj with no k passing qc
  SummExpDR.obj.no.pass <- run_multi_k_functions(x = test_obj,
                                                analysis_name = 'clust.analysis',
                                                redDimKey = NULL,
                                                dims_use = NULL,
                                                k_use = k_use,
                                                num_iter = num_iter,
                                                algorithm = algorithm,
                                                pct_resample = pct_resample,
                                                affi_K = affi_K,
                                                sigma = sigma,
                                                base_seed = base_seed,
                                                thresh_pval = 0.1,
                                                thresh_ARI = 1.0,
                                                thresh_proportion = 1.0)
  expect_error(get_soln_k(object = SummExpDR.obj.no.pass,
                          pick_best_valid = T,
                          k = NULL,
                          add_2_metadata = F,
                          analysis_use = 1))

})
