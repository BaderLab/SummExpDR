###############################################################################
context('test DNAm_anno')
devtools::load_all()
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# illumina annotations of EPIC probes
epic_anno <- minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
# epic_anno$CHR <- paste0('chr', epic_anno$CHR)
set.seed(12345L)
epic_anno <- epic_anno[sample(1:nrow(epic_anno), size = 25),]
# setup transcript ranges for checks
transcripts <- GenomicFeatures::transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)
names(transcripts) <- transcripts$tx_name
tx_name <- names(transcripts)
seq_names <- as.vector(GenomicRanges::seqnames(transcripts), mode = 'character')
neg_strand <- as.vector(GenomicRanges::strand(transcripts), mode = 'character') == '-'
starts <- as.vector(GenomicRanges::start(transcripts), mode = 'integer')
ends <- as.vector(GenomicRanges::end(transcripts), mode = 'integer')
inds_use <- 1:nrow(epic_anno)
# setup granges_prep list to test annotate_probes function
promoter_prep <- create_granges_prep(transcripts, operation = 'upstream', analysis_name = 'tss200', bp = 200L)
gene_body_prep <- create_granges_prep(transcripts, operation = 'in', analysis_name = 'gene_body')
granges_prep_list <- list(promoter_prep, gene_body_prep)

# setup points to use for testing
chr1_pos_idx <- sample(which(!neg_strand & seq_names == 'chr1'), size = 1, replace = FALSE) # index of gene on positive strand, chr1
chr1_neg_idx <- sample(which(neg_strand & seq_names == 'chr1'), size = 1, replace = FALSE) # index of gene on negative strand, chr1
tryCatch({
  stopifnot(length(chr1_pos_idx) == 1)
  stopifnot(length(chr1_neg_idx) == 1)
}, error = function(e) {
  stop('failed to find index for chr1 positive strand transcript and or chr1 negative strand transcript')
})
start_end1 <- c(starts[chr1_pos_idx], ends[chr1_pos_idx])
start_end2 <- c(starts[chr1_neg_idx], ends[chr1_neg_idx])
upstream_1 <- start_end1[1] - 150 # positive strand case, point is within 200bp upstream of tss
out_1a <- start_end1[1] - 201 # positive strand case, too far upstream of tss (> 200bp)
out_1b <- start_end1[2] + 1 # positive strand case, outside of gene (further than transcript end point)
within_1 <- start_end1[1] + 10 # after TSS, before end
upstream_2 <- start_end2[2] + 150 # negative strand case, point is within 200bp upstream of tss
out_2a <- start_end2[2] + 201 # negative strand case, point is too far upstream of tss
out_2b <- start_end2[2] - 1 # negative strand case, point is beyond TSS
# declare global variables for run options for annotate_probes
chr_col <- 'chr'
coord_col <- 'pos'

# Run Tests

test_that('.check_upstream', {
  # check on positive and negative strands,
  expect_identical(TRUE, .check_upstream(upstream_1, chr = 'chr1', regions = transcripts, bp = 200)[chr1_pos_idx])
  expect_identical(FALSE, .check_upstream(out_1a, chr = 'chr1', regions = transcripts, bp = 200)[chr1_pos_idx])
  expect_identical(FALSE, .check_upstream(out_1b, chr = 'chr1', regions = transcripts, bp = 200)[chr1_pos_idx])
  expect_identical(FALSE, .check_upstream(within_1, chr = 'chr1', regions = transcripts, bp = 200)[chr1_pos_idx])
  expect_identical(TRUE, .check_upstream(upstream_2, chr = 'chr1', regions = transcripts, bp = 200)[chr1_neg_idx])
  expect_identical(FALSE, .check_upstream(out_2a, chr = 'chr1', regions = transcripts, bp = 200)[chr1_neg_idx])
  expect_identical(FALSE, .check_upstream(out_2b, chr = 'chr1', regions = transcripts, bp = 200)[chr1_neg_idx])

  # test results on whole set of transcripts
  upstream_pos <- ((starts - upstream_1) <= 200) & ((starts - upstream_1) > 0) & (seq_names == 'chr1')
  upstream_neg <- ((upstream_1 - ends) <= 200) & ((upstream_1 - ends) > 0) & (seq_names == 'chr1')
  check_upstream_result <- .check_upstream(upstream_1, chr = 'chr1', regions = transcripts, bp = 200)
  expect_identical(upstream_pos[!neg_strand], check_upstream_result[!neg_strand])
  expect_identical(upstream_neg[neg_strand], check_upstream_result[neg_strand])
  # anything not on chr1 should not be labeled as upstream
  expect_true(all(check_upstream_result[seq_names != 'chr1'] == FALSE))

})

test_that('.check_isin', {
  # do for individual coordinates for 1 gene
  expect_identical(TRUE, .check_isin(within_1, chr = 'chr1', regions = transcripts)[chr1_pos_idx])
  expect_identical(FALSE, .check_isin(out_1b, chr = 'chr1', regions = transcripts)[chr1_pos_idx])
  expect_identical(FALSE, .check_isin(upstream_1, chr = 'chr1', regions = transcripts)[chr1_pos_idx])
  within_idx <- (within_1 >= starts & within_1 <= ends) & seq_names == 'chr1'
  # do for whole set of ranges
  within_results <- .check_isin(within_1, chr = 'chr1', regions = transcripts)
  expect_identical(within_idx, within_results)
  expect_true(all(within_results[seq_names != 'chr1'] == FALSE))
})

test_that('annotate_probes: serial', {
  # run the annotation function
  epic_anno_new <- annotate_probes(illumina_anno = epic_anno,
                                   granges_prep_list = granges_prep_list,
                                   prefix = 'su2c_gencode_v12',
                                   coord_column = coord_col,
                                   chr_column = chr_col,
                                   n_cores = 1)
  # check gene bodies
  reference_genes <- vapply(inds_use, FUN.VALUE = character(1), FUN = function(x) {
    chr_x <- epic_anno[x, chr_col]
    coord_x <- epic_anno[x, coord_col]
    mapped_gene_inds <- .check_isin(coord = coord_x, chr = chr_x, regions = transcripts)
    mapped_genes <- names(transcripts)[mapped_gene_inds]
    if (length(mapped_genes) == 0) {
      mapped_genes <- ''
    } else {
      mapped_genes <- paste(mapped_genes, collapse = ';')
    }
    return(mapped_genes)
  })
  gene_body_annot <- epic_anno_new$su2c_gencode_v12_gene_body
  expect_identical(reference_genes, gene_body_annot)
  # check promoters
  reference_promoters <- vapply(inds_use, FUN.VALUE = character(1), FUN = function(x) {
    chr_x <- epic_anno[x, chr_col]
    coord_x <- epic_anno[x, coord_col]
    mapped_gene_inds <- .check_upstream(coord = coord_x, chr = chr_x, regions = transcripts, bp = 200)
    mapped_genes <- names(transcripts)[mapped_gene_inds]
    if (length(mapped_genes) == 0) {
      mapped_genes <- ''
    } else {
      mapped_genes <- paste(mapped_genes, collapse = ';')
    }
    return(mapped_genes)
  })
  promoter_annot <- epic_anno_new$su2c_gencode_v12_tss200
  expect_identical(reference_promoters, promoter_annot)
})

test_that('annotate_probes: parallel functionality', {
  # run the annotation function
  epic_anno_serial <- annotate_probes(illumina_anno = epic_anno,
                                      granges_prep_list = granges_prep_list,
                                      prefix = 'su2c_gencode_v12',
                                      coord_column = coord_col,
                                      chr_column = chr_col,
                                      n_cores = 1)
  # run annotation function, parallel
  epic_anno_parallel <- annotate_probes(illumina_anno = epic_anno,
                                        granges_prep_list = granges_prep_list,
                                        prefix = 'su2c_gencode_v12',
                                        coord_column = coord_col,
                                        chr_column = chr_col,
                                        n_cores = 3)
  expect_identical(epic_anno_serial, epic_anno_parallel)
})

devtools::unload('GenomicFeatures')
devtools::unload('TxDb.Hsapiens.UCSC.hg19.knownGene')
