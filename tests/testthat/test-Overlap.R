context('test overlap utils')

test_that('fisher_test_sets produces sane output', {
  set.a <- letters[1:10]
  set.b <- letters[2:12]
  result.1 <- fisher_test_sets(set1 = set.a, set2 = set.b, master_set = letters)

  expect_lt(result.1$p.value, 0.01)

  set.seed(20190128)
  set.a.random <- sample(letters, size = 10)
  set.b.random <- sample(letters, size = 11)
  result.2 <- fisher_test_sets(set1 = set.a.random, set2 = set.b.random, master_set = letters)

  expect_gt(result.2$p.value, 0.10)


})

test_that('do_fisher_multi produces expected output', {
  abc <- letters[1:3]
  bcd <- letters[2:4]
  def <- letters[4:6]
  abcd <- letters[1:4]
  efgh <- letters[5:8]
  hijk <- letters[8:11]
  xyza <- c('x', 'y', 'z', 'a')
  test.list1.orig <- list(abc = abc, bcd = bcd, def = def,
                          abcd = abcd, efgh = efgh, hijk = hijk, xyza = xyza)

  query.orig <- c('r', 'q', 'a', 'd', 'c')

  for (rnd_seed in 42:44) {
    set.seed(rnd_seed)
    test.list1 <- test.list1.orig[sample(1:length(test.list1.orig), length(test.list1.orig), replace = FALSE)]
    query <- query.orig[sample(1:length(query.orig), length(query.orig), replace = FALSE)]
    fisher_multi_result <- do_fisher_multi(query, test.list1, master_set = letters, FDR = 1.0, n_cores = 1)
    # check that FDR correction done properly
    expect_equal(fisher_multi_result$FDR, p.adjust(fisher_multi_result$p_val, method = 'BH'))
    # check that p-values and intersection calculated correctly
    for (pathway_name in names(test.list1)) {
      fisher_result <- fisher_test_sets(query, test.list1[[pathway_name]], master_set = letters)
      fisher_multi_row <- fisher_multi_result[fisher_multi_result$pathway_name == pathway_name,]
      # check p-value calculated correctly
      expect_equal(fisher_multi_row$p_val, fisher_result$p.value)
      # check intersection made correctly
      expect_equal(sort(unlist(strsplit(fisher_multi_row$intersect, split = ';'))), sort(intersect(query, test.list1[[pathway_name]])))
      # check that pathway contents stored
      expect_equal(unlist(strsplit(fisher_multi_row$pathway, split = ';')), test.list1[[pathway_name]])
    }
  }
})

test_that('do_fisher_multi produces expected output in parallel case', {
  abc <- letters[1:3]
  bcd <- letters[2:4]
  def <- letters[4:6]
  abcd <- letters[1:4]
  efgh <- letters[5:8]
  hijk <- letters[8:11]
  xyza <- c('x', 'y', 'z', 'a')
  test.list1.orig <- list(abc = abc, bcd = bcd, def = def,
                          abcd = abcd, efgh = efgh, hijk = hijk, xyza = xyza)

  query.orig <- c('r', 'q', 'a', 'd', 'c')

  for (rnd_seed in 42:44) {
    set.seed(rnd_seed)
    test.list1 <- test.list1.orig[sample(1:length(test.list1.orig), length(test.list1.orig), replace = FALSE)]
    query <- query.orig[sample(1:length(query.orig), length(query.orig), replace = FALSE)]
    fisher_multi_result <- do_fisher_multi(query, test.list1, master_set = letters, FDR = 1.0, n_cores = 1)
    # check that p-values and intersection calculated correctly
    for (n_cores in 2:3) {
      fisher_multi_result_parallel <- do_fisher_multi(query, test.list1, master_set = letters, FDR = 1.0, n_cores = n_cores)
      expect_identical(fisher_multi_result, fisher_multi_result_parallel)
    }
  }
})

test_that('calc_overlap produces expected output', {
  abc <- letters[1:3]
  bcd <- letters[2:4]
  def <- letters[4:6]

  test.list1.orig <- list(abc = abc, bcd = bcd, def = def)

  abcd <- letters[1:4]
  efgh <- letters[5:8]
  hijk <- letters[8:11]
  xyza <- c('x', 'y', 'z', 'a')

  test.list2.orig <- list(abcd = abcd, efgh = efgh, hijk = hijk, xyza = xyza)
  for (rnd_seed in 42:44) {
    # previous version of calc_overlap produced incorrect output if list names were not ordered.
    # permuting lists
    set.seed(rnd_seed)
    test.list1 <- test.list1.orig[sample(1:length(test.list1.orig), length(test.list1.orig), replace = FALSE)]
    test.list2 <- test.list2.orig[sample(1:length(test.list2.orig), length(test.list2.orig), replace = FALSE)]
    for (metric in c('num.ovr', 'jaccard', 'ovr.coef', 'fisher.p')) {
      if (metric == 'fisher.p') {
        master_set <- letters
      } else {
        master_set <- NULL
      }
      M1 <- calc_overlap(list1 = test.list1,
                         list2 = NULL,
                         metric = metric,
                         n_cores = 1,
                         master_set = master_set)
      M2 <- calc_overlap(list1 = test.list1,
                         list2 = test.list2,
                         metric = metric,
                         n_cores = 1,
                         master_set = master_set)
      ## ref.1a = abc x bcd
      ## ref.1b = bcd x def
      ## ref.2a = abc x abcd
      ## ref.2b = def x efgh
      if (metric == 'num.ovr') {
        ref.1a <- 2
        ref.1b <- 1
        ref.2a <- 3
        ref.2b <- 2
      } else if (metric == 'jaccard') {
        ref.1a <- 1/2
        ref.1b <- 1/5
        ref.2a <- 3/4
        ref.2b <- 2/5
      } else if (metric == 'ovr.coef') {
        ref.1a <- 2/3
        ref.1b <- 1/3
        ref.2a <- 1
        ref.2b <- 2/3
      } else if (metric == 'fisher.p') {
        ref.1a <- -log10(fisher_test_sets(set1 = abc,
                                          set2 = bcd,
                                          master_set = letters)$p.value)
        ref.1b <- -log10(fisher_test_sets(set1 = bcd,
                                          set2 = def,
                                          master_set = letters)$p.value)
        ref.2a <- -log10(fisher_test_sets(set1 = abc,
                                          set2 = abcd,
                                          master_set = letters)$p.value)
        ref.2b <- -log10(fisher_test_sets(set1 = def,
                                          set2 = efgh,
                                          master_set = letters)$p.value)
      }
      expect_equal(M1['abc','bcd'], ref.1a)
      expect_equal(M1['bcd', 'def'], ref.1b)
      expect_equal(M2['abc', 'abcd'], ref.2a)
      expect_equal(M2['def', 'efgh'], ref.2b)

      ## expect M1 is mxm, expect M2 is mxn
      expect_equal(dim(M1), c(length(test.list1), length(test.list1)))
      expect_equal(dim(M2), c(length(test.list1), length(test.list2)))
    }
  }


})

test_that('calc_overlap parallel functionality works', {
  # baderlab_human_pathways_dec_2018 is a list of character vectors
  gene2probe <- baderlab_human_pathways_dec_2018[grep('NFKB|TGF|FGF|EGF|WNT', names(baderlab_human_pathways_dec_2018))]
  gene2probe_a <- gene2probe[1:40]
  gene2probe_b <- gene2probe[60:101]
  for (metric in c('num.ovr', 'jaccard', 'ovr.coef', 'fisher.p')) {
    if (metric == 'fisher.p') {
      master_set <- unique(unlist(gene2probe))
    } else {
      master_set <- NULL
    }
    expect_equal(calc_overlap(gene2probe_a, gene2probe_b, n_cores = 1, metric = metric, master_set = master_set),
                 calc_overlap(gene2probe_a, gene2probe_b, n_cores = 2, metric = metric, master_set = master_set))
    expect_equal(calc_overlap(gene2probe_a, gene2probe_b, n_cores = 1, metric = metric, master_set = master_set),
                 calc_overlap(gene2probe_a, gene2probe_b, n_cores = 3, metric = metric, master_set = master_set))
    # case of longer list2, is switching of inputs and transposing of output handled correctly
    expect_equal(calc_overlap(gene2probe_b, gene2probe_a, n_cores = 1, metric = metric, master_set = master_set),
                 calc_overlap(gene2probe_b, gene2probe_a, n_cores = 2, metric = metric, master_set = master_set))
  }
})
