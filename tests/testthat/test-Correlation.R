set.seed(42L)
# 100 x 100 matrix
data1 <- matrix(c(rnorm(9000), rnorm(1000, 3)), nrow = 100)
rownames(data1) <- paste0('feat', 1:nrow(data1))
colnames(data1) <- paste0('sample_', as.character(1:ncol(data1)))
data1 <- t(data1)
# 100 x 100 matrix, with randomly selected values from same distributions
data2 <- matrix(c(rnorm(9000), rnorm(1000, 3)), nrow = 100)
rownames(data2) <- paste0('feat', 1:nrow(data1))
colnames(data2) <- paste0('sample_', as.character(1:ncol(data1)))
data2 <- t(data2)

testthat::test_that('get_cor_vals works', {
  serial_result <- get_cor_vals(data1, data2, n_cores = 1L, fdr_filter = 1.0)
  testthat::expect_true(nrow(serial_result) == nrow(data1)*nrow(data2))
  serial_result_self_only <- get_cor_vals(data1, data2, n_cores = 1L, self_only = TRUE, fdr_filter = 1.0)
  testthat::expect_true(all(serial_result_self_only$feat_x == serial_result_self_only$feat_y))
  ref_diff <- cor(data1[,'feat1'], data2[,'feat2'])
  ref_self <- cor(data1[,'feat5'], data2[,'feat5'])
  res_diff <- serial_result[serial_result$feat_x == 'feat1' & serial_result$feat_y == 'feat2', ]
  res_self <- serial_result_self_only[serial_result_self_only$feat_x == 'feat5', ]
  # test correlation values calculated
  testthat::expect_equal(ref_diff, res_diff$cor)
  testthat::expect_equal(ref_self, res_self$cor)
  # repeat for parallel
  parallel_result <- get_cor_vals(data1, data2, n_cores = 2L, fdr_filter = 1.0)
  testthat::expect_equal(serial_result, parallel_result)
  parallel_result_self_only <- get_cor_vals(data1, data2, n_cores = 2L, self_only = TRUE, fdr_filter = 1.0)
  testthat::expect_equal(serial_result_self_only, parallel_result_self_only)
})
