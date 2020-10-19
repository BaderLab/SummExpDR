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

test_that('calc_overlap produces expected output', {
  abc <- letters[1:3]
  bcd <- letters[2:4]
  def <- letters[4:6]

  test.list1 <- list(abc = abc, bcd = bcd, def = def)

  abcd <- letters[1:4]
  efgh <- letters[5:8]
  hijk <- letters[8:11]
  xyza <- c('x', 'y', 'z', 'a')

  test.list2 <- list(abcd = abcd, efgh = efgh, hijk = hijk, xyza = xyza)
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
