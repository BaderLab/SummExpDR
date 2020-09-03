##########################################################################

#' Run GSVA
#' @include SummExpDR.R MultiExpIntegration.R
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @param x
#' @param assay assay to pull data from
#' @param genesets named list of genesets to pass to GSVA::gsva
#' @param expt_use experiment to use if this is a multiExp object
#' @value return a SummarizedExperiment, SummExpDR, or multiExp object with colData updated with gsva scores

setGeneric('run_gsva', function(x, assay, genesets, expt_use = NULL, ...) standardGeneric('run_gsva'))

setMethod('run_gsva',
          signature = 'SummarizedExperiment',
          function(x, assay, genesets, expt_use = NULL, ...) {
            exp_data <- SummarizedExperiment::assay(x, assay)
            gsva_scores <- GSVA::gsva(exp_data, gset.idx.list = genesets, ...)
            SummarizedExperiment::colData(x) <- cbind(SummarizedExperiment::colData(x),
                                                      S4Vectors::DataFrame(t(gsva_scores)[colnames(x), ]))
            return(x)
          })

setMethod('run_gsva',
          signature = 'SummExpDR',
          function(x, assay, genesets, expt_use = NULL, ...) {
            summ_exp <- getSummExp(x)
            x@summ_exp <- run_gsva(summ_exp, assay, genesets, expt_use = NULL, ...)
            return(x)
          })

setMethod('run_gsva',
          signature = 'multiExp',
          function(x, assay, genesets, expt_use = NULL, ...) {
            tryCatch({stopifnot(is.character(expt_use) && length(expt_use == 1))},
                     error = function(e) {
                       stop('expt_use must be character of length 1 for signature multiExp')
                     })
            # note that feature ids or mapped back to original name in input data matrix
            summ_exp <- getSummExp(x)
            row_data <- SummarizedExperiment::rowData(summ_exp)
            row_data_expt <- row_data[row_data$expt == expt_use, ]
            name_mapping <- row_data_expt$orig_id
            names(name_mapping) <- row_data_expt$feat_id
            summ_exp_subs <- summ_exp[names(name_mapping), ]
            rownames(summ_exp_subs) <- name_mapping
            exp_data <- SummarizedExperiment::assay(summ_exp_subs, assay)
            gsva_scores <- GSVA::gsva(exp_data, gset.idx.list = genesets, ...)
            rownames(gsva_scores) <- paste(rownames(gsva_scores), expt, sep = '_')
            SummarizedExperiment::colData(summ_exp) <- cbind(SummarizedExperiment::colData(summ_exp),
                                                      S4Vectors::DataFrame(t(gsva_scores)[colnames(summ_exp), ]))
            x@summ_exp <- summ_exp
            return(x)
          })
