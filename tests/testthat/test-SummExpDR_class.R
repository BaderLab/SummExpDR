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
  se <- SummarizedExperiment(assays = list(assay1 = data_mat), rowData = row_data, colData = col_data)
  SEDR <- create_SummExpDR(se)
  testthat::expect_error(create_SummExpDR(data_mat))
})


data_mat2 <- matrix(rnorm(2000), nrow = 10)
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
                                                              size = ncol(data_mat2),
                                                              replace = TRUE),
                                             id = sample(7:9,
                                                         size = ncol(data_mat2),
                                                         replace = TRUE),
                                             row.names = colnames(data_mat2),
                                             stringsAsFactors = FALSE))

testthat::test_that('runPCA works', {
  se2 <- SummarizedExperiment(assays = list(assay1 = data_mat2), rowData = row_data2, colData = col_data2)
  SEDR2 <- create_SummExpDR(se2)
  SEDR2 <- runPCA(SEDR2, 1, '_assay1')
  loadings_mat <- getLoadings(SEDR2)
})
