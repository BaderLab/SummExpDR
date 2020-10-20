###############################################################################
#############################################
### RNA differential expression functions ###

#' Run differential expression analysis
#'
#' @param expt_data SummarizedExperiment object with RNA-seq data or microarray
#' data, or SummExpDR with said SummarizedExperiment object in @summ_exp slot. If using
#' limma as selected differential expression pipeline, can handle DNAm m-values as well
#' or any data matrix where features are expected to have a Gaussian distribution
#' @param class_type_use class denoted in metadata (must be one of metadata columns)
#' @param class1 value for first class, which will be 'numerator' for log fold change
#' @param class2 will be 'denominator' for log fold change. if set to NULL,
#' defined as the complement of class1 (i.e. all other classes), and called 'all.others'
#' @param covariates to use in linear model or GLM if using 'voom', 'limma', or
#' 'DESEQ2' pipelines. Must be character. If length 0, design matrix consists only
#' of intercept and covariate of interest (class_type_use).
#' If length 1 or more, design matrix consists if intercept, covariate of interest,
#' ass less as additional covariates specified by this parameter
#' @param pipeline one of voom, limma, or DESEQ2.
#' @param assay_use which assay to pull out of SummarizedExperiment object.
#' should be 1 unless there are multiple assays and you want to get an assay that isn't
#' the first 1
#' @param output_dir top output directory for
#' @param prefix specified prefix for rank files (i.e. filename = [prefix].rnk). If null,
#' defaults to [class1]_vs_[class2]
#' @param overwrite overwrite output files of prior run
#' @param write_output write output to file, i.e. diff exp result, rankfile. files will be written
#' in subdirectory of output_dir
#' @param ... other arguments to pass to run_limma or run_DESEQ2
#' @return
#' If input is SummarizedExperiment, list containing following elements:
#' 1) diff_exp = results of differential expression analysis using voom +
#' limma  if 'voom' is argument for pipeline(see \code{\link{run_limma}}),
#' limma if 'limma' selected for pipeline (see \code{\link{run_limma}}) or
#' DESEQ2 (see {\code{\link{run_DESEQ2}}) if 'DESEQ2' provided.
#' 2) rnk = ranked genelist ranked by LFC*(-log10(adjusted p.val)). if using 'voom' pipeline,
#' LFC = log2FC
#' 3) comparison = character vector, length 1, with string formatted as [class_type_use]_[class1]_vs_[class2].
#' note that reported LFC between class1 and class2 is always positive for class1 upregulated genes
#' and negative for class2 downregulated genes
#' If input is su2c_obj class object, return su2c_obj with above list with results for
#' input's SummarizedExperiment object assigned as a named entry to the diff.phen attribute
#' @export
#'
#' @examples
#'

run_diff_exp <- function(expt_data,
                         class_type_use,
                         class1,
                         class2 = NULL,
                         covariates = character(0),
                         assay_use = 1,
                         pipeline = 'voom',
                         write_output = T,
                         output_dir = './',
                         prefix = '',
                         overwrite = F,
                         ...) {
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  if (class(expt_data) == 'SummExpDR') {
    obj <- expt_data
    expt_data <- getSummExp(obj)
    de_result <- run_diff_exp(expt_data = expt_data,
                              class_type_use = class_type_use,
                              class1 = class1,
                              class2 = class2,
                              covariates = covariates,
                              assay_use = assay_use,
                              pipeline = pipeline,
                              write_output = F,
                              prefix = prefix,
                              output_dir = output_dir,
                              ...)
    name.result <- de_result$comparison
    obj <- setAnalyses(obj, name.result, de_result)
    if (write_output) {
      write_RNA_DE_output(obj = obj,
                          result_name = name.result,
                          top_output_dir = output_dir,
                          overwrite = overwrite)
    }
    return(obj)
  } else if (class(expt_data) == 'SummarizedExperiment') {
    stopifnot(is.character(prefix) && (length(prefix) == 1))
    comparison.string <- paste0(prefix, class_type_use, '_', class1, '_v_')
    if (is.null(class2)) {
      comparison.string <- paste0(comparison.string, 'AllOthers')
    } else {
      comparison.string <- paste0(comparison.string, class2)
    }
    if (pipeline %in% c('voom', 'limma')) {
      if (pipeline == 'voom') {
        use_voom <- T
      } else {
        use_voom <- F
      }
      output<- run_limma(expt_data = expt_data,
                         assay_use = assay_use,
                         class_type_use = class_type_use,
                         class1 = class1,
                         class2 = class2,
                         covariates = covariates,
                         use_voom = use_voom,
                         ...)

      LFC.padj <- output[,c('logFC', 'adj.P.Val')]
      ranked.genelist <- rank_and_order(LFC_padj_df = LFC.padj,
                                        write.rankfile = F)

    } else if (pipeline == 'DESEQ2') {
      output <- run_DESEQ2(expt_data = expt_data,
                           class_type_use = class_type_use,
                           class1 = class1,
                           class2 = class2,
                           assay_use = assay_use,
                           covariates = covariates,
                           ...)
      LFC.padj <- output[,c('log2FoldChange', 'padj')]
      ranked.genelist <- rank_and_order(LFC_padj_df = LFC.padj,
                                        write.rankfile = F)
    } else {
      stop(paste('pipeline', pipeline, 'not implemented for this function.',
                 'valid pipelines include voom, limma, and DESEQ2'))
    }
    list.return <- list(diff_exp = output,
                        rnk = ranked.genelist,
                        comparison = comparison.string)
    if (write_output) {
      write_RNA_DE_output(obj = list.return,
                          top_output_dir = output_dir,
                          overwrite = overwrite)
    }
    return(list.return)
  } else {
    stop(paste('function not implemented for class', class(expt_data)))
  }
}

#' Write ranked genelist and diff exp results to output files
#'
#' @param obj either a su2c_obj class object or a list output by run_diff_exp
#' @param result_name name of result to retrieve from obj@diff.phen. NULL by default,
#' must be specified if obj is of class su2c_obj
#' @param top_output_dir top output directory below which a subdirectory containing
#' analysis output is kept
#' @param overwrite overwrite contents of output directory (subdirectory of top_output_dir).
#' In practice,  a system call is made to remove output directory and its contents
#' before making it again if output directory exists. If set to FALSE, will raise error if
#' directory already exists. If set to TRUE, and in the unlikely event the subdirectory created has
#' the same name as something you don't want to remove, provided that the user running the
#' process has the appropriate permissions, the directory will be removed.
#' @details if list output by run_diff_exp, write out ranked genelist entry (rnk),
#' and write out differential expression entry (diff_exp)
#' @return Nothing returned, just write out ranked genelist (see rank and order) and
#' differential expression results (see run_voom).
#' @export
#'
#' @examples
write_RNA_DE_output <- function(obj, result_name = NULL, top_output_dir = './', overwrite = F) {
  if (class(obj) == 'SummExpDR') {
    result <- getAnalyses(obj, result_name)
    write_RNA_DE_output(obj = result, top_output_dir = top_output_dir, overwrite = overwrite)
  } else if (class(obj) == 'list') {
    if (!all(names(obj) == c('diff_exp', 'rnk', 'comparison'))) {
      stop('unexpected list input. should be output of run_diff_exp for SummarizedExperiment class')
    }
    stopifnot(dir.exists(top_output_dir))

    diff_exp <- obj$diff_exp
    rnk <- obj$rnk
    comparison <- obj$comparison

    ## BE EXTREMELY CAREFUL ABOUT BELOW LINE AND IF STATEMENT.
    ## AVOID DOING SOMETHING LIKE system(paste('rm -rf', top_output_dir))
    subdir <- file.path(top_output_dir, comparison)
    if (dir.exists(subdir)) {
      if (overwrite) {
        warning(paste('Directory', subdir, 'exists,
                      removing since overwrite set to TRUE'))
        system(paste('rm -rf', subdir))
      } else {
        stop(paste('Directory', subdir, 'exists, overwrite set to FALSE'))
      }
    }
    if (!dir.exists(subdir)) {
      dir.create(subdir)
    }
    write.table(diff_exp, file = file.path(subdir, paste0(comparison, '_diff_exp.tsv')),
                sep = '\t', row.names = T, quote = F)
    write.table(rnk, file = file.path(subdir, paste0(comparison, '.rnk')),
                sep = '\t', row.names = F, quote = F)
  }
}
###
#' Subset sample metadata by class
#' @inheritParams run_diff_exp
#' @param metadata metadata
#' @param class_type_use
#' @param class1
#' @param class2
#' @return metadata subsetted by class if 2 classes specified, or metadata
#' with column class_type_use changed to have classes [class1] and 'all.others'
#' @export
#'
#' @examples
subs_by_class <- function(metadata,
                          class_type_use,
                          class1,
                          class2) {
  if (is(metadata, 'DataFrame')) {
    metadata <- as.data.frame(metadata)
  }
  stopifnot(is.data.frame(metadata))
  metadata[, class_type_use] <- as.character(metadata[, class_type_use])
  if (class_type_use %in% colnames(metadata)) {
    meta.col <- metadata[, class_type_use]
  } else {
    stop(paste('metadata has no feature', class_type_use))
  }
  if (!class1 %in% meta.col) {
    stop('class1 must be in metadata column')
  }
  if (is.null(class2)) {
    class2 <- 'all.others'
    metadata[meta.col != class1, class_type_use] <- class2
  } else  {
    class2 <- as.character(class2)
    if (class2 %in% meta.col) {
      metadata <- metadata[meta.col %in% c(class1, class2),]
    }
    else {
      stop(paste('class2 must be null or must be in metadata column', class_type_use))
    }
  }
  ## for any columns with factors, make sure that they are updated to remove
  ## any categories which no longer exist in the subsetted metadata
  for (i in 1:ncol(metadata)) {
    if (is.factor(metadata[,i])) {
      metadata[,i] <- factor(as.character(metadata[,i]))
    }
  }
  ## convert class of interest to factor
  metadata[,class_type_use] <- factor(metadata[,class_type_use], levels = c(class2, class1))
  return(metadata)
}

#' Run limma or voom and limma
#' @inheritParams run_diff_exp
#' @param expt_data
#' @param class_type_use
#' @param class1
#' @param class2
#' @param assay_use assay to use from expt_data. defaults to first slot.
#' @param covariates covariates to include in linear model
#' @param use_voom use edgeR/voom/limma workflow if TRUE. If FALSE, use limma differential
#' expression experimental data provided. Should set to TRUE for RNA-seq data,
#' FALSE for log RMA values for microarray data
#' @param remove_low_cts
#' @return differential expression results table as output by
#' toptable
#' @export
#'
#' @examples
run_limma <- function(expt_data,
                      class_type_use,
                      class1,
                      class2,
                      assay_use = 1,
                      covariates = character(0),
                      use_voom = F,
                      remove_low_cts = T) {
  ## get data
  data.mat <- SummarizedExperiment::assay(expt_data, assay_use)

  meta.data <- SummarizedExperiment::colData(expt_data)
  ## subset according to comparison performed
  meta.data <- subs_by_class(metadata = meta.data,
                             class_type_use = class_type_use,
                             class1 = class1,
                             class2 = class2)
  data.mat <- data.mat[,rownames(meta.data)]
  ## make design matrix
  formula.string <- paste('~', class_type_use)
  stopifnot(is.character(covariates))
  if (length(covariates) >= 1) {
    for (i in 1:length(covariates)) {
      cov.i <- covariates[i]
      if (! cov.i %in% colnames(meta.data)) {
        stop(paste('covariate', cov.i, 'not found in sample metadata'))
      } else if (cov.i == class_type_use) {
        warning(paste(cov.i, 'in specified covariates, but is already variable of interest'))
      } else {
        formula.string <- paste(formula.string, '+', cov.i)
      }
    }
  }
  formula.use <- formula(formula.string)
  mod.frame <- model.frame(formula.use, data = meta.data)
  design.mat <- model.matrix(formula.use, data = mod.frame)
  if (use_voom) {
    ## voom + limma used for RNA-seq data
    if (!all(is.integer(data.mat))) {
      stop('expect count data from expt_data for assay chosen')
    }
    ## make DGElist and do TMM normalization
    DGE <- edgeR::DGEList(counts = data.mat)
    if (remove_low_cts) {
      keep <- edgeR::filterByExpr(DGE)
    }
    counts <- DGE$counts
    counts <- counts[keep,]
    counts <- counts[,rownames(meta.data)]
    DGE <- edgeR::DGEList(counts = counts)
    DGE <- edgeR::calcNormFactors(DGE)
    ## run diff exp testing with voom + limma
    lm.input <- limma::voom(DGE, design.mat, plot = F)
  } else {
    lm.input <- data.mat
  }

  ## fit linear model, run diff exp
  fit <- limma::lmFit(lm.input, design.mat)
  fit <- limma::eBayes(fit)
  coef.name <- paste0(class_type_use, class1)
  if (!coef.name %in% colnames(fit$coefficients)) {
    stop(paste('expected name of coefficient for class1 to be ',
               coef.name, 'in model returned'))
  }
  results <- limma::topTable(fit,
                             coef = coef.name,
                             number = Inf,
                             sort.by = 'logFC')
  return(results)
}

#' Run DESEQ2 on count data
#' @inheritParams run_diff_exp
#' @param expt_data SummarizedExperiment object with count data.
#' @param assay_use index or name of assay in expt_data to pass to DESeq2. Should
#' indicate a matrix of counts, and an error will be raised by \code{\link{DESeq2::DESeqDataSet}}
#' if the matrix accessed is not a matrix of counts.
#' @param class_type_use
#' @param class1
#' @param class2
#' @param n.cores number of cores to use for DESeq2 differential expression
#' @param assay_use which assay to use from expt_data. should indicate counts
#' @param covariates
#' @param remove_low_cts remove low counts? If this is the case, remove all genes
#' with row sums less than rowsum.thresh
#' @param rowsum.thresh see remove_low_cts
#' @return data frame of results with LFC defined by log2(c1/c2), where c1 and
#' c2 give the expression of class1 and class2, respectively. Uses BH corrected
#' p-values, 2 sided test (altHypothesis = 'greaterAbs'). Uses default
#' arguments values for DESeq2::\code{\link{results}}
#' @export
#'
#' @examples
run_DESEQ2 <- function(expt_data,
                       class_type_use,
                       class1,
                       class2 = NULL,
                       n.cores = 1,
                       assay_use = 1,
                       covariates = character(0),
                       remove_low_cts = F,
                       rowsum.thresh = 5) {
  if (!is(expt_data, 'SummarizedExperiment')) {
    stop('expect expt_data to be SummarizedExperiment')
  }

  if (n.cores == 0) {
    n.cores <- parallel::detectCores()
  }
  if (n.cores == 1) {
    use.parallel <- F
  } else if (n.cores > 1) {
    use.parallel <- T
  } else {
    stop('n.cores should be >= 1 or specified as 0 for automatic detection')
  }
  if (remove_low_cts) {
    ## by default assumed that low counts have been removed previously,
    ## but including as a convenience option
    low.cts.ind <- which(apply(SummarizedExperiment::assay(expt_data, assay_use),
                               MARGIN = 1,
                               FUN = sum) < rowsum.thresh)
    expt_data <- expt_data[-low.cts.ind,]
  }

  meta.data <- SummarizedExperiment::colData(expt_data)
  class.labels <- as.character(meta.data[,class_type_use])
  if (!class1 %in% class.labels) {
    stop(paste('argument class1:', class1, 'not found in metadata column', class_type_use, 'for
               colData of expt_data'))
  }
  if (is.null(class2)) {
    class2 <- 'all.others'
    class.labels[class.labels != class1] <- class2
  } else if (!class2 %in% class.labels) {
    stop(paste('argument class2:', class2, 'not found in metadata column', class_type_use, 'for
               colData of expt_data. specify class2 as NULL or valid class for this category'))
  }

  meta.data[,class_type_use] <- factor(class.labels)
  ## get data
  data.mat <- SummarizedExperiment::assay(expt_data, assay_use)
  expt_data <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = data.mat),
                                                          colData = meta.data)


  formula.string <- paste('~', class_type_use)
  stopifnot(is.character(covariates))
  if (length(covariates) >= 1) {
    for (i in 1:length(covariates)) {
      cov.i <- covariates[i]
      if (! cov.i %in% colnames(meta.data)) {
        stop(paste('covariate', cov.i, 'not found in sample metadata'))
      } else if (cov.i == class_type_use) {
        warning(paste(cov.i, 'in specified covariates, but is already variable of interest'))
      } else {
        formula.string <- paste(formula.string, '+', cov.i)
      }
    }
  }
  dds <- DESeq2::DESeqDataSet(se = expt_data, design = formula(formula.string))
  dds <- DESeq2::estimateSizeFactors(dds)
  dds <- DESeq2::DESeq(dds, parallel = use.parallel, BPPARAM = BiocParallel::MulticoreParam(n.cores, exportglobals = FALSE))
  DESEQ2.result <- DESeq2::results(dds,
                                   contrast = c(class_type_use, class1, class2),
                                   altHypothesis = 'greaterAbs',
                                   pAdjustMethod = 'BH')
  DESEQ2.result <- filter_NA(DESEQ2.result)
  return(DESEQ2.result)
}

#' Filter NA values from deseq2 result
#'
#' @param DESEQ2.result DESeqResults object
#'
#' @return DESeqResults object with NA values removed
#'
#'
#' @examples
filter_NA <- function(DESEQ2.result) {
  if (is(DESEQ2.result, 'DESeqResults')) {
    DESEQ2.result <- DESEQ2.result[which(!is.na(DESEQ2.result$padj)),]
  } else {
    stop('expect object of class DESeqResults')
  }
}


### Rank scoring functions

#' Rank and order genes according to Log Fold Change and adjusted p-value
#'
#' @param LFC_padj_df data.frame with log fold change RNA exp in column 1, FDR (padj) in second
#' @param write.rankfile write a rankfile?
#' @param output_dir filepath to directory where you want to write rankfile
#' @param fname filename, to be appended to out_
#' @details scores are calculated as -log10(padj)*sign(LFC). To account for Inf and -Inf values,
#' any scores with abs(score) == Inf are converted to sign(LFC)*(max(abs(score.i)) + runif(0,1)),
#' s.t. score.i != inf. runif(0,1) random uniformly distributed variable, from range 0,1
#' @return data frame with genes in first column, score in second
#' @export
#'
#' @examples
rank_and_order <- function(LFC_padj_df,
                           write.rankfile = T,
                           output_dir = getwd(),
                           fname = 'ranked_list.rnk') {
  # rank each gene in provided dataframe by the product of logfold change and -log10(p_value)
  rank_score <- function(LFC, p.adj) {
    return(sign(LFC)*(-log10(p.adj)))
  }
  rnk_scores <- apply(LFC_padj_df, MARGIN = 1, FUN = function(x) {
    LFC <- x[1]
    p.adj <- x[2]
    return(rank_score(LFC, p.adj))
  })
  # Take into account log10(0) = -Inf
  # get largest magninude Non-Inf rnk_score
  rnk_scores_max_abs <- max(abs(rnk_scores[ !is.infinite(abs(rnk_scores))]))
  set.seed(12345)
  # assign max value*sign + unif dist for any infinite elements
  rnk_scores_non_inf <- sapply( rnk_scores,
                                function(x) replace(x, is.infinite(abs(x)),
                                                    sign(x) * (rnk_scores_max_abs + runif(n = 1))) )

  LFC_padj_df$rank_score <- rnk_scores_non_inf
  LFC_padj_rnk_df <- LFC_padj_df[order(rnk_scores_non_inf, decreasing = T),]
  rnk_df <- data.frame(genes = rownames(LFC_padj_rnk_df),
                       rank_score = LFC_padj_rnk_df$rank_score,
                       stringsAsFactors = F)
  if(write.rankfile) {
    write.table(rnk_df, file = file.path(output_dir, fname),
                sep = '\t',
                quote = F, row.names = F)
  }
  return(rnk_df)
}
