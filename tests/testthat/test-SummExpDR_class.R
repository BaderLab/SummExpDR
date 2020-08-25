testthat::context('SummExpRD_class')

data_mat <- matrix(rnorm(9), nrow = 3)
rownames(data_mat) <- c('a', 'b', 'c')
colnames(data_mat) <- c('x', 'y', 'z')
col_data <- S4Vectors::DataFrame(data.frame(color = c('red', 'blue', 'green'),
                                            int = 1:3, row.names = colnames(data_mat),
                                            stringsAsFactors = FALSE))
row_data <- S4Vectors::DataFrame(data.frame(texture = c('soft', 'medium', 'rough'),
                                            id = 7:9, row.names = rownames(data_mat),
                                            stringsAsFactors = FALSE))

testthat::test_that('class created correctly', {
  se <- SummarizedExperiment::SummarizedExperiment(assays = list(assay1 = data_mat), rowData = row_data, colData = col_data)
  SEDR <- create_SummExpDR(se)
  testthat::expect_error(create_SummExpDR(data_mat))
})


data_mat2 <- t(matrix(rnorm(2000), nrow = 100, ncol = 20))
rownames(data_mat2) <- paste0('feat', 1:nrow(data_mat2))
colnames(data_mat2) <- paste0('sample', 1:ncol(data_mat2))
col_data2 <- S4Vectors::DataFrame(data.frame(color = sample(c('red', 'blue', 'green'),
                                                            size = ncol(data_mat2),
                                                            replace = TRUE),
                                             int = sample(1:3,
                                                          size = ncol(data_mat2),
                                                          replace = TRUE),
                                             row.names = colnames(data_mat2),
                                             stringsAsFactors = FALSE))
row_data2 <- S4Vectors::DataFrame(data.frame(texture = sample(c('soft', 'medium', 'rough'),
                                                              size = nrow(data_mat2),
                                                              replace = TRUE),
                                             id = sample(7:9,
                                                         size = nrow(data_mat2),
                                                         replace = TRUE),
                                             row.names = rownames(data_mat2),
                                             stringsAsFactors = FALSE))
se2 <- SummarizedExperiment::SummarizedExperiment(assays = list(assay1 = data_mat2), rowData = row_data2, colData = col_data2)
SEDR2 <- create_SummExpDR(se2)
SEDR2 <- runPCA(SEDR2, 'assay1', '_assay1')
pca_ref <- prcomp(t(SummarizedExperiment::assay(se2, 'assay1')), scale. = TRUE, center = TRUE)

testthat::test_that('runPCA works', {
  reduced_dims <- getReducedDims(SEDR2, key = 'PCA_assay1')
  loadings_mat <- getLoadings(reduced_dims)
  coords_mat <- getEmbeddings(reduced_dims)
  testthat::expect_equal(dim(loadings_mat), c(20,20))
  testthat::expect_equal(dim(coords_mat), c(20, 100))
  # check calculated coordinates

  testthat::expect_true(all(t(pca_ref$x) == coords_mat))
})

testthat::test_that('varianceExplained works', {
  # in this case calculating variance explained by first 3 PCs
  dims_use <- 1:3
  var_expl_res <- varianceExplained(SEDR2, key = 'PCA_assay1', dims_use = dims_use)
  total_var_ref <- sum(pca_ref$sdev^2)
  r2_dim_ref <- (pca_ref$sdev^2)[dims_use]/total_var_ref
  r2_total_ref <- sum((pca_ref$sdev^2)[dims_use])/total_var_ref
  names(r2_dim_ref) <- paste0('PC', dims_use)
  testthat::expect_equal(var_expl_res$r2_by_dim, r2_dim_ref)
  testthat::expect_equal(var_expl_res$r2_all, r2_total_ref)
})

testthat::test_that('runVarimax works', {
  SEDR2 <- runVarimax(SEDR2, key = 'PCA_assay1', dims_use = 1:3)
  vmax_result <- getReducedDims(SEDR2, 'varimax')
  vmax_result_loadings <- getLoadings(vmax_result)
  pca_result <- getReducedDims(SEDR2, 'PCA_assay1')
  pca_loadings <- getLoadings(pca_result, red_dims = 1:3)
  vmax_ref <- varimax(t(pca_loadings))
  vmax_ref_loadings <- t(vmax_ref$loadings[,])
  rownames(vmax_ref_loadings) <- paste0('VM', 1:3)
  testthat::expect_equal(vmax_ref_loadings, vmax_result_loadings)
})
