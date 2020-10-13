###############################################################################
context('DNAm_utils')
# set up variables for all tests
set.seed(42)
# setup smaller gene probe mapping for speed
all_map_probes <- unlist(strsplit(illumina_EPIC_hg19_10b4_gene_2_probe_tss200$mapped_probes, split = ';'))
common_probes <- intersect(rownames(DNAm_bval), all_map_probes)
diff_probes <- setdiff(all_map_probes, rownames(DNAm_bval))
gene_2_probe <- illumina_EPIC_hg19_10b4_gene_2_probe_tss200
probes_use <- paste(c(sample(common_probes, 100), sample(diff_probes, 400)), collapse = '|')
gene_2_probe <- gene_2_probe[grep(probes_use, gene_2_probe$mapped_probes), ]
# make tempdir for tests
tmpdir <- './DNAm_tests_tmp/'
testthat::setup({
  
  if (dir.exists(tmpdir)) {
    system(paste('rm -rf', tmpdir))
  }
  dir.create(tmpdir)
})


# run tests
test_that('calculate_probe_avg works', {

  probe_avg_result <- calculate_probe_avg(DNAm_bval, gene_2_probe, n_cores = 1)
  DNAm_bval_probe_avg <- probe_avg_result[[1]]
  probes_per_gene <- probe_avg_result[[2]]
  expect_equal(names(probes_per_gene), rownames(DNAm_bval_probe_avg))

  has_probe_in_data <- check_min_probes(1, data_mat = DNAm_bval, gene_probe_mapping = gene_2_probe)
  genes_w_probes <- gene_2_probe$gene[has_probe_in_data]
  genes_w_no_probes <- gene_2_probe$gene[-which(has_probe_in_data)]
  # make sure that cases of genes with + without probes are tested
  stopifnot(length(genes_w_probes > 1) && length(genes_w_no_probes > 1))
  for (i in 1:10) {
    # test for genes w probes. mean bval of probes should be same
    gene_test <- sample(genes_w_probes, 1)
    probes_use <- get_mapped_probes(gene_test, gene_2_probe)
    probes_use <- intersect(probes_use, rownames(DNAm_bval))
    if (length(probes_use) == 1) {
      avg_bval_test <- DNAm_bval[probes_use,]
    } else {
      avg_bval_test <- Matrix::colMeans(DNAm_bval[probes_use,])
    }
    expect_equal(avg_bval_test, DNAm_bval_probe_avg[gene_test, ])
  }
  for (i in 1:10) {
    # test for genes w/o probes. these genes should nbot be in rows of final output
    gene_test <- sample(genes_w_no_probes, 1)
    expect_false(gene_test %in% rownames(DNAm_bval_probe_avg))
  }
})

test_that('calculate_probe_avg_works, parallel', {
  probe_avg_res_serial <- calculate_probe_avg(DNAm_bval, gene_2_probe, n_cores = 1)
  probe_avg_res_parallel <- calculate_probe_avg(DNAm_bval, gene_2_probe, n_cores = 2)
  expect_identical(probe_avg_res_serial, probe_avg_res_parallel)
})


rm(gene_2_probe)

testthat::teardown({
  if (dir.exists(tmpdir)) {
    system(paste('rm -rf', tmpdir))
  }
})
