testthat::context('MultiExpIntegration')
set.seed(42L)
# 100 x 100 matrix
data1 <- matrix(c(rnorm(9000), rnorm(1000, 3)), nrow = 100)
rownames(data1) <- paste0('data1_', 1:nrow(data1))
colnames(data1) <- as.character(1:ncol(data1))
# 60 x 200 matrix
data2 <- matrix(c(rnorm(2400, 5, 2), rnorm(9600, 2, 5)), nrow = 60)
rownames(data2) <- paste0('data2_', 1:nrow(data2))
colnames(data2) <- as.character(51:250)
# 50 samples (columns) shared between matrices

colData1 <- S4Vectors::DataFrame(SampleID = colnames(data1))
colData2 <- S4Vectors::DataFrame(SampleID = colnames(data2))

summExp1 <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = data1), colData = colData1)
summExp2 <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = data2), colData = colData2)


testthat::test_that('multiExp Class Created Correctly', {
  multiExp <- createMultiExp(summ_exp_list = list(data1 = summExp1, data2 = summExp2), assays_use = c(data1 = 'mat', data2 = 'mat'))
  summ_exp <- getSummExp(multiExp)
  assay_data <- SummarizedExperiment::assay(summ_exp, 'stacked')
  testthat::expect_equal(dim(assay_data), c(160, 250))
  num_na <- sum(is.na(assay_data))
  # 60 features x 50 samples with data1 but no data2
  # 100 features x 150 samples with data2 but no data1
  expect_na <- 60*50 + 100*150
  testthat::expect_equal(num_na, expect_na)
  row_data <- SummarizedExperiment::rowData(summ_exp)
  testthat::expect_equal(sort(row_data[row_data$expt == 'data1', 'orig_id']), sort(rownames(data1)))
  testthat::expect_equal(sort(row_data[row_data$expt == 'data2', 'orig_id']), sort(rownames(data2)))
})

testthat::test_that('multiExp data imputation produces expected output', {
  multiExp <- createMultiExp(summ_exp_list = list(data1 = summExp1, data2 = summExp2), assays_use = c(data1 = 'mat', data2 = 'mat'))
  multiExp <- imputeExpData(multiExp)
  summ_exp <- getSummExp(multiExp)
  expect_na <- 60*50 + 100*150
  imputed_data <- SummarizedExperiment::assay(summ_exp, 'imputed_mat')
  testthat::expect_equal(sum(is.na(imputed_data)), 0)
  is_imputed <- SummarizedExperiment::assay(summ_exp, 'is_imputed')
  testthat::expect_equal(sum(is_imputed), expect_na)
  testthat::expect_lt(abs(mean(imputed_data[101:160,as.character(1:50)]) - 5), 1)
})

testthat::test_that('data scaling works properly', {
  multiExp <- createMultiExp(summ_exp_list = list(data1 = summExp1, data2 = summExp2), assays_use = c(data1 = 'mat', data2 = 'mat'))
  multiExp <- imputeExpData(multiExp)
  multiExp <- scaleExpData(multiExp, assay = 'imputed_mat')
  summ_exp <- getSummExp(multiExp)
  scaled_data <- SummarizedExperiment::assay(summ_exp, 'scaled')
  testthat::expect_true(all(matrixStats::rowMeans2(scaled_data) - 0 < .Machine$double.eps))
  covar_mat <- (scaled_data %*% t(scaled_data))/(ncol(scaled_data) - 1)
  var_data1 <- sum(diag(covar_mat)[1:100])
  var_data2 <- sum(diag(covar_mat)[101:160])
  testthat::expect_equal(var_data1, 1L)
  testthat::expect_equal(var_data2, 1L)
})
#
# testthat::test_that('runPCA workflow works', {
#
# })
