.onLoad <- function(libname, pkgname) {
  tryCatch({
    SKLEARN_IMPUTE <<- reticulate::import('sklearn.impute')
    SummarizedExperiment <<- getFromNamespace('SummarizedExperiment', 'SummarizedExperiment')
  }, error = function(e) {
    stop('Failed to load sklearn.impute via reticulate. check if python3, sklearn installed.')
  })
}
