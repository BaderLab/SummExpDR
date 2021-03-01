context('DE_utils')
data_dir <- '../../data'
load(file.path(data_dir, 'RNA_SE.RData'))
load(file.path(data_dir, 'cavalli_2017_testdata.RData'))

tmpdir <- './tmp_RNA_DE'

testthat::setup({
  dir.create(tmpdir)
})



## Construct Test Data --------------------------------------------------------

# For microarray data, we will be adding a simmulated batch effect
# and a corresponding batch variable
set.seed(54321)
cavalli.data <- SummarizedExperiment::assay(cavalli_2017_testdata, 1)
samples.mutate <- sample(colnames(cavalli.data), floor(ncol(cavalli.data)/2),
                         replace = F)
num.genes <- nrow(cavalli.data)
batch.effect <- matrix(rnorm(n = num.genes*length(samples.mutate), mean = 3, sd = 0.5),
                       nrow = num.genes, ncol = length(samples.mutate))
cavalli.data[,samples.mutate] <- cavalli.data[,samples.mutate] + batch.effect
SummarizedExperiment::assay(cavalli_2017_testdata, 1) <- cavalli.data
cavalli.batch <- rep(1, ncol(cavalli.data))
cavalli.batch[colnames(cavalli.data) %in% samples.mutate] <- 2
cavalli_2017_testdata$batch <- cavalli.batch
rm(cavalli.data, samples.mutate, num.genes, batch.effect, cavalli.batch)
cavalli_2017_testdata$subgroup <- cavalli_2017_testdata$subgroup.ch1

# Add Another set of classes to RNA-seq data
RNA_SE$KC_Sex <- paste(RNA_SE$Sex, as.character(RNA_SE$KC), sep = '_')

## metadata subsetteing -------------------------------------------------------
context('metadata_subsetting DE_utils')
test_that('metadata subsetting works', {
  RNA_data <- cavalli_2017_testdata
  ## case of 1 vs all
  meta_data <-  SummarizedExperiment::colData(RNA_data)
  result_1 <- subs_by_class(metadata = meta_data,
                            class_type_use = 'subgroup',
                            class1 = 'WNT',
                            class2 = NULL)
  result_1_wnt <- result_1[result_1$subgroup == 'WNT', ]
  expect_equal(as.character(result_1_wnt[, 'subgroup']),
               as.character(meta_data[rownames(result_1_wnt), 'subgroup']))
  expect_true(all(result_1[result_1$subgroup != 'WNT', 'subgroup'] == 'all.others'))
  expect_equal(rownames(result_1), rownames(meta_data))
  ## case of class 1 vs class 2
  result_2 <- subs_by_class(metadata = meta_data,
                            class_type_use = 'subgroup',
                            class1 = 'WNT',
                            class2 = 'Group3')
  result_2_wnt <- result_2[result_2$subgroup == 'WNT', ]
  result_2_non_wnt <- result_2[result_2$subgroup != 'WNT', ]
  expect_equal(as.character(result_2_wnt[, 'subgroup']),
               as.character(meta_data[rownames(result_2_wnt), 'subgroup']))
  expect_true(all(as.character(result_2_non_wnt[,'subgroup']) == 'Group3'))
  expect_equal(as.character(result_2_non_wnt[, 'subgroup']),
               as.character(meta_data[rownames(result_2_non_wnt), 'subgroup']))
  expect_equal(nrow(result_2), sum(as.numeric(meta_data$subgroup %in% c('WNT', 'Group3'))))
  ## check that when subsetting metadata, all factors columns are updated
  ## to remove non-existent categories
  meta_data_subs <- meta_data[meta_data$subgroup %in% c('WNT', 'Group3', 'Group4'),]
  result_3 <- subs_by_class(metadata = meta_data_subs,
                            class_type_use = 'subgroup',
                            class1 = 'Group3',
                            class2 = 'Group4')
  for (i in 1:ncol(result_3)) {
    if (is.factor(result_3[,i])) {
      expect_true(all(levels(result_3[,i]) %in% result_3[,i]))
    }
  }
})
## run limma ------------------------------------------------------------------
context('run_limma')
test_that('run_limma  voom mode', {
  RNA_data <- RNA_SE

  ## Test with 1 vs all comparison
  ## no covariates included
  result <- run_limma(expt_data = RNA_data,
                      class_type_use = 'KC_Sex',
                      class1 = 'M_TRUE',
                      class2 = NULL,
                      covariates = character(0),
                      assay_use = 1,
                      use_voom = T,
                      remove_low_cts = T)
  ## covariates included
  result.cov <- run_limma(expt_data = RNA_data,
                          class_type_use = 'KC_Sex',
                          class1 = 'M_TRUE',
                          class2 = NULL,
                          covariates = c('Age'),
                          assay_use = 1,
                          use_voom = T,
                          remove_low_cts = T)
  # expect inclusion of covariates to change result
  expect_false(identical(rownames(result), rownames(result.cov)))
  ## Try with only 2 classes
  result.2.class <- run_limma(expt_data = RNA_data,
                          class_type_use = 'KC_Sex',
                          class1 = 'M_TRUE',
                          class2 = 'M_FALSE',
                          covariates = c('Age'),
                          assay_use = 1,
                          use_voom = T,
                          remove_low_cts = T)
  # expect doing class1 v class2 comparison causes different result
  expect_false(identical(rownames(result.cov), rownames(result.2.class)))
  ## make sure direction of LFC is correct
  result.KC <- run_limma(expt_data = RNA_data,
                        class_type_use = 'KC',
                        class1 = TRUE,
                        class2 = FALSE,
                        covariates = c('Age'),
                        assay_use = 1,
                        use_voom = T,
                        remove_low_cts = T)
  KC_TRUE <- colnames(RNA_data)[RNA_data$KC == TRUE]
  KC_FALSE <- colnames(RNA_data)[RNA_data$KC == FALSE]
  RNA_mat <- SummarizedExperiment::assay(RNA_data, 1)
  result.KC <- result.KC[result.KC$adj.P.Val < 0.05, ]
  class1_genes <- rownames(result.KC)[result.KC$logFC > 0]
  class2_genes <- rownames(result.KC)[result.KC$logFC < 0]
  # take mean of DE genes for each class, compare
  expect_gt(mean(RNA_mat[class1_genes, KC_TRUE]), mean(RNA_mat[class1_genes, KC_FALSE]))
  expect_lt(mean(RNA_mat[class2_genes, KC_TRUE]), mean(RNA_mat[class2_genes, KC_FALSE]))
})

test_that('run_limma works under default (microarray) mode', {
  RNA_data <- cavalli_2017_testdata
  data.mat <- SummarizedExperiment::assay(RNA_data, 1)
  batch1.samples <- colnames(RNA_data)[RNA_data$batch == 1]
  non.batch1.samples <- colnames(RNA_data)[RNA_data$batch == 2]
  wnt.samples <- colnames(RNA_data)[RNA_data$subgroup == 'WNT']
  nonwnt.samples <- colnames(RNA_data)[RNA_data$subgroup != 'WNT']
  gr3.samples <- colnames(RNA_data)[RNA_data$subgroup == 'Group3']
  ## Test with 1 vs all comparison
  ## no covariates included
  result <- run_limma(expt_data = RNA_data,
                      class_type_use = 'subgroup',
                      class1 = 'WNT',
                      class2 = NULL,
                      covariates = character(0),
                      assay_use = 1,
                      use_voom = F)
  ## covariates included
  result.cov <- run_limma(expt_data = RNA_data,
                          class_type_use = 'subgroup',
                          class1 = 'WNT',
                          class2 = NULL,
                          covariates = 'batch',
                          assay_use = 1,
                          use_voom = F)
  # expect inclusion of covariates to change result
  expect_false(identical(rownames(result), rownames(result.cov)))
  ## no covariates included, batch 1 samples only
  result.batch1 <- run_limma(expt_data = RNA_data[,batch1.samples],
                             class_type_use = 'subgroup',
                             class1 = 'WNT',
                             class2 = NULL,
                             covariates = character(0),
                             assay_use = 1,
                             use_voom = F)
  ## batch 1 samples only, vs group 3
  results.gr3 <- run_limma(expt_data = RNA_data[,batch1.samples],
                           class_type_use = 'subgroup',
                           class1 = 'WNT',
                           class2 = 'Group3',
                           covariates = character(0),
                           assay_use = 1,
                           use_voom = F)
  # expect doing class1 v class2 comparison causes different result
  expect_false(identical(rownames(result.batch1), rownames(results.gr3)))
  ## make sure direction of LFC is correct
  upreg.genes.cov <- rownames(result.cov[result.cov$logFC > 0
                                         & result.cov$adj.P.Val < 0.05,])
  upreg.genes.batch1 <- rownames(result.batch1[result.batch1$logFC > 0
                                               & result.batch1$adj.P.Val < 0.05,])
  upreg.genes.v.gr3 <- rownames(results.gr3[results.gr3$logFC > 0
                                            & results.gr3$adj.P.Val < 0.05,])
  expect_gt(length(intersect(upreg.genes.batch1, upreg.genes.cov)), 0)
  wnt.batch1 <- intersect(wnt.samples, batch1.samples)
  gr3.batch1 <- intersect(gr3.samples, batch1.samples)
  nonwnt.batch1 <- intersect(nonwnt.samples, batch1.samples)
  # Wnt upregulated genes are higher in WNT when batch held constant
  expect_gt(mean(data.mat[upreg.genes.batch1, wnt.batch1]),
            mean(data.mat[upreg.genes.batch1, nonwnt.batch1]))
  expect_gt(mean(data.mat[upreg.genes.cov, wnt.batch1]),
            mean(data.mat[upreg.genes.cov, nonwnt.batch1]))
  # Wnt upregulated genes vs group 3 higher in Wnt than in group 3
  expect_gt(mean(data.mat[upreg.genes.v.gr3, wnt.batch1]),
            mean(data.mat[upreg.genes.v.gr3, gr3.batch1]))
})
## run_DESEQ2 -----------------------------------------------------------------
context('DESEQ2')

test_that('filter_NA produces correct output, fails on incorrect input', {
  RNA_data <- RNA_SE
  RNA_data$KC <- factor(RNA_data$KC)
  ## make DESEQ2 results
  dds <- DESeq2::DESeqDataSet(RNA_data, design = formula('~ KC'))
  dds <- DESeq2::estimateSizeFactors(dds)
  res <- DESeq2::results(DESeq2::DESeq(dds), contrast = c('KC', 'TRUE', 'FALSE'),
                         altHypothesis = 'greaterAbs', pAdjustMethod = 'BH')
  ## filter
  # stop if no NAs as test becomes pointless
  stopifnot(any(is.na(res$padj)))
  res.filt <- filter_NA(res)
  expect_equal(res.filt, res[which(!is.na(res$padj)),])
  expect_true(all(is.numeric(res.filt$padj)))
  expect_true(all(rownames(res.filt) %in% rownames(res)))

  ## test for failures
  expect_error(filter_NA('asdf'))
  expect_error(filter_NA(as.data.frame(res.filt)))
  expect_error(filter_NA(213))
})

test_that('run_DESEQ2 produces expected output with and without covariates', {
  RNA_data <- RNA_SE
  RNA_data$KC <- factor(RNA_data$KC)

  ## test without covariates
  dds.no.cov <- DESeq2::DESeqDataSet(RNA_data, design = formula('~ KC'))
  dds.no.cov <- DESeq2::estimateSizeFactors(dds.no.cov)
  ref.no.cov <- DESeq2::results(DESeq2::DESeq(dds.no.cov), contrast = c('KC', 'TRUE', 'FALSE'),
                                altHypothesis = 'greaterAbs', pAdjustMethod = 'BH')
  ref.no.cov <- filter_NA(ref.no.cov)

  test.no.cov <- run_DESEQ2(expt_data = RNA_data,
                            class_type_use = 'KC',
                            class1 = 'TRUE',
                            class2 = 'FALSE',
                            assay_use = 1,
                            covariates = character(0),
                            remove_low_cts = F,
                            n_cores = 1)

  expect_identical(test.no.cov, ref.no.cov)

  ## test with covariates
  dds.cov <- DESeq2::DESeqDataSet(RNA_data, design = formula('~ KC + Age + Sex'))
  dds.cov <- DESeq2::estimateSizeFactors(dds.cov)
  ref.cov <- DESeq2::results(DESeq2::DESeq(dds.cov), contrast = c('KC', 'TRUE', 'FALSE'),
                             altHypothesis = 'greaterAbs', pAdjustMethod = 'BH')
  ref.cov <- filter_NA(ref.cov)

  test.cov <- run_DESEQ2(expt_data = RNA_data,
                         class_type_use = 'KC',
                         class1 = 'TRUE',
                         class2 = 'FALSE',
                         assay_use = 1,
                         covariates = c('Age', 'Sex'),
                         remove_low_cts = F,
                         n_cores = 1)

  expect_identical(test.cov, ref.cov)
  expect_false(identical(test.cov, test.no.cov))
})

test_that('run_DESEQ2 handles class2 arguments correctly', {
  RNA_data <- RNA_SE
  ref.metadata <- data.frame(SummarizedExperiment::colData(RNA_data))
  RNA.ref.data <- SummarizedExperiment::SummarizedExperiment(assays = SummarizedExperiment::assays(RNA_data),
                                                             colData = ref.metadata)
  RNA.ref.data$KC_Sex[RNA_data$KC_Sex != 'M_TRUE'] <- 'all.others'
  ## test without class2 specified
  dds.null <- DESeq2::DESeqDataSet(RNA.ref.data, design = formula('~ KC_Sex + Age'))
  dds.null <- DESeq2::estimateSizeFactors(dds.null)
  ref.null <- DESeq2::results(DESeq2::DESeq(dds.null), contrast = c('KC_Sex', 'M_TRUE', 'all.others'),
                              altHypothesis = 'greaterAbs', pAdjustMethod = 'BH')
  ref.null <- filter_NA(ref.null)

  test.null <- run_DESEQ2(expt_data = RNA_data,
                          class_type_use = 'KC_Sex',
                          class1 = 'M_TRUE',
                          class2 = NULL,
                          assay_use = 1,
                          covariates = 'Age',
                          remove_low_cts = F,
                          n_cores = 1)

  expect_identical(test.null, ref.null)

  ## test with class2 specified
  dds.class2 <- DESeq2::DESeqDataSet(RNA_data, design = formula('~ KC_Sex + Age'))
  dds.class2 <- DESeq2::estimateSizeFactors(dds.class2)
  ref.class2 <- DESeq2::results(DESeq2::DESeq(dds.class2), contrast = c('KC_Sex', 'M_TRUE', 'M_FALSE'),
                                altHypothesis = 'greaterAbs', pAdjustMethod = 'BH')
  ref.class2 <- filter_NA(ref.class2)

  test.class2 <- run_DESEQ2(expt_data = RNA_data,
                            class_type_use = 'KC_Sex',
                            class1 = 'M_TRUE',
                            class2 = 'M_FALSE',
                            assay_use = 1,
                            covariates = 'Age',
                            remove_low_cts = F,
                            n_cores = 1)

  expect_identical(test.class2, ref.class2)
  expect_false(identical(test.class2, test.null))
})
## rank and order -------------------------------------------------------------
context('rank and order')
test_that('rank_and_order produces expected output', {
  RNA_data <- RNA_SE
  result <- run_limma(expt_data = RNA_data,
                      class_type_use = 'KC',
                      class1 = 'TRUE',
                      class2 = NULL,
                      assay_use = 1,
                      use_voom = T,
                      remove_low_cts = T)
  LFC.padj <- result[,c('logFC', 'adj.P.Val')]
  ## manually setting some adjusted p.values to 0 to give case
  ## of infinite adjsuted p values
  set.seed(12345)
  LFC.padj$adj.P.Val[sample(1:nrow(LFC.padj), size = 10)] <- 0

  rnk1 <- rank_and_order(LFC_padj_df = LFC.padj,
                         write.rankfile = F,
                         output_dir = tmpdir,
                         fname = 'result1.rnk')
  ## make sure no file output if write.rankfile = F
  expect_true(!any(grepl('result1.rnk', dir(tmpdir))))
  ## testing contents of rnk1
  lfc1 <- LFC.padj[rnk1$genes, 'logFC']
  p.adj <- LFC.padj[rnk1$genes, 'adj.P.Val']
  raw.score <- sign(lfc1)*(-log10(p.adj))
  inf.vals <- which(abs(raw.score) == Inf)
  raw.score.non.inf <- raw.score[-inf.vals]
  expect_equal(raw.score.non.inf, rnk1$rank_score[-inf.vals])
  if (length(inf.vals) > 0) {
    rnk1.inf <- rnk1[inf.vals, 'rank_score']
    lfc1.inf <- lfc1[inf.vals]
    max.val <- max( abs(raw.score[-inf.vals]) )
    ## absolute values of rank score should be bounded
    ## between max absolute value for raw scores and that value + 1
    expect_true( all( abs(rnk1.inf) >= rep(max.val , length(rnk1.inf) ) ) )
    expect_true( all( abs(rnk1.inf) <= rep( max.val + 1,  length(rnk1.inf) ) ) )
    ## sign of score should be same as that of corresponding LFC
    expect_true( all(sign(rnk1.inf) == sign(lfc1[inf.vals])))
  } else {
    stop('expected at least some infinite rank scores in manually calculated scores')
  }
  rnk2 <- rank_and_order(LFC_padj_df = LFC.padj,
                         write.rankfile = T,
                         output_dir = tmpdir,
                         fname = 'result2.rnk')
  ## changing write.rankfile should not change output
  expect_identical(rnk1, rnk2)
  ## rankfile contents should match output R dataframe
  rnk2.load <- read.table(file = file.path(tmpdir, 'result2.rnk'),
                          sep = '\t',
                          header = T,
                          stringsAsFactors = F)
  expect_equal(rnk2, rnk2.load)
})

## run_diff_exp ---------------------------------------------------------------
context('run_diff_exp')

test_that('run_diff_exp, DESEQ2', {

  # note that creating SummExpDR object results in renaming according to R string naming conventions, which can mess with equivalency tests
  # if testing for identical rownames, so we take out SummarizedExperiment object once that's done
  SummExpDR_obj <- create_SummExpDR(RNA_SE)
  RNA_data <- getSummExp(SummExpDR_obj)
  class1 <- 'M_TRUE'
  assay_use <- 1

  for (class2 in c('M_FALSE', 'NULL')) {
    if (class2 == 'NULL') {
      class2 <- NULL
    }
    for (covar in list(none = character(0), Age = 'Age')) {

      DESEQ2.result <- run_DESEQ2(expt_data = RNA_data,
                                  class_type_use = 'KC_Sex',
                                  class1 = class1,
                                  class2 = class2,
                                  assay_use = 1,
                                  covariates = covar)

      DE.result.SE <- run_diff_exp(expt_data = RNA_data,
                                   class_type_use = 'KC_Sex',
                                   class1 = class1,
                                   class2 = class2,
                                   covariates = covar,
                                   assay_use = 1,
                                   pipeline = 'DESEQ2',
                                   write_output = F,
                                   output_dir = './',
                                   prefix = '',
                                   overwrite = F)


      DE.result.SummExpDR.obj <- run_diff_exp(expt_data = SummExpDR_obj,
                                               class_type_use = 'KC_Sex',
                                               class1 = class1,
                                               class2 = class2,
                                               covariates = covar,
                                               assay_use = 1,
                                               pipeline = 'DESEQ2',
                                               write_output = F,
                                               output_dir = './',
                                               prefix = '',
                                               overwrite = F)

      expect_equal(DE.result.SE$diff_exp, DESEQ2.result)
      expect_equal(DE.result.SE$rnk, rank_and_order(DESEQ2.result[,c('log2FoldChange', 'padj')] ,
                                                    write.rankfile = F))
      analysis_key <- getAnalyses_keys(DE.result.SummExpDR.obj)
      expect_identical(getAnalyses(DE.result.SummExpDR.obj, key = analysis_key), DE.result.SE)
      expect_true(all(is.finite(DE.result.SE$rnk$rank_score)))
    }
  }

})

test_that('run_diff_exp, voom', {


  class1 <- 'M_TRUE'
  assay_use <- 1
  pipeline <- 'voom'
  # note that creating SummExpDR object results in renaming according to R string naming conventions, which can mess with equivalency tests
  # if testing for identical rownames, so we take out SummarizedExperiment object once that's done
  SummExpDR_obj <- create_SummExpDR(RNA_SE)
  RNA_data <- getSummExp(SummExpDR_obj)

  for (class2 in c('M_FALSE', 'NULL')) {
    if (class2 == 'NULL') {
      class2 <- NULL
    }
    for (covar in list(none = character(0), Age = 'Age')) {

      limma.result <- run_limma(expt_data = RNA_data,
                                  class_type_use = 'KC_Sex',
                                  class1 = class1,
                                  class2 = class2,
                                  assay_use = 1,
                                  covariates = covar,
                                  use_voom = TRUE)

      DE.result.SE <- run_diff_exp(expt_data = RNA_data,
                                   class_type_use = 'KC_Sex',
                                   class1 = class1,
                                   class2 = class2,
                                   covariates = covar,
                                   assay_use = 1,
                                   pipeline = pipeline,
                                   write_output = F,
                                   output_dir = './',
                                   prefix = '',
                                   overwrite = F)


      DE.result.SummExpDR.obj <- run_diff_exp(expt_data = SummExpDR_obj,
                                              class_type_use = 'KC_Sex',
                                              class1 = class1,
                                              class2 = class2,
                                              covariates = covar,
                                              assay_use = 1,
                                              pipeline = pipeline,
                                              write_output = F,
                                              output_dir = './',
                                              prefix = '',
                                              overwrite = F)

      expect_equal(DE.result.SE$diff_exp, limma.result)
      expect_equal(DE.result.SE$rnk, rank_and_order(limma.result[,c('logFC', 'adj.P.Val')] ,
                                                    write.rankfile = F))
      analysis_key <- getAnalyses_keys(DE.result.SummExpDR.obj)
      expect_identical(getAnalyses(DE.result.SummExpDR.obj, key = analysis_key), DE.result.SE)
      expect_true(all(is.finite(DE.result.SE$rnk$rank_score)))
    }
  }

})

test_that('run_diff_exp, limma', {

  class1 <- 'WNT'
  assay_use <- 1
  pipeline <- 'limma'
  # note that creating SummExpDR object results in renaming according to R string naming conventions, which can mess with equivalency tests
  # if testing for identical rownames, so we take out SummarizedExperiment object once that's done
  SummExpDR_obj <- create_SummExpDR(cavalli_2017_testdata)
  RNA_data <- getSummExp(SummExpDR_obj)

  for (class2 in c('Group3', 'NULL')) {
    if (class2 == 'NULL') {
      class2 <- NULL
    }
    for (covar in list(none = character(0), batch = 'batch')) {

      limma.result <- run_limma(expt_data = RNA_data,
                                class_type_use = 'subgroup',
                                class1 = class1,
                                class2 = class2,
                                covariates = covar,
                                use_voom = F)


      DE.result.SE <- run_diff_exp(expt_data = RNA_data,
                                   class_type_use = 'subgroup',
                                   class1 = class1,
                                   class2 = class2,
                                   covariates = covar,
                                   assay_use = 1,
                                   pipeline = pipeline,
                                   write_output = F,
                                   output_dir = './',
                                   prefix = '',
                                   overwrite = F)

      DE.result.SummExpDR.obj <- run_diff_exp(expt_data = SummExpDR_obj,
                                             class_type_use = 'subgroup',
                                             class1 = class1,
                                             class2 = class2,
                                             covariates = covar,
                                             assay_use = 1,
                                             pipeline = pipeline,
                                             write_output = F,
                                             output_dir = './',
                                             prefix = '',
                                             overwrite = F)

      expect_equal(DE.result.SE$diff_exp, limma.result)
      expect_equal(DE.result.SE$rnk, rank_and_order(limma.result[,c('logFC', 'adj.P.Val')] ,
                                                    write.rankfile = F))
      analysis_key <- getAnalyses_keys(DE.result.SummExpDR.obj)
      expect_identical(getAnalyses(DE.result.SummExpDR.obj, key = analysis_key), DE.result.SE)
      expect_true(all(is.finite(DE.result.SE$rnk$rank_score)))
    }
  }
})


testthat::teardown({
  system(paste('rm -rf', tmpdir))
})
