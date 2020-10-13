set.seed(42L)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(SummExpDR)
library(GenomicFeatures)
library(minfi)
# assemble small DNAm b-val dataset
DNAm_data <- readRDS('~/projects/su2c_v2/data/preprocessed/HGCC/DNAm_SE.rds')
DNAm_bval <- SummarizedExperiment::assay(DNAm_data, 'DNAm_bval')
DNAm_bval <- DNAm_bval[-which(apply(DNAm_bval, MARGIN = 1, FUN = function(x) {any(is.na(x))})), ]
DNAm_bval <- DNAm_bval[sample(rownames(DNAm_bval), 5000), ]
save(DNAm_bval, file = 'data/DNAm_bval.RData')

# assemble gene to probe mapping for genes to probes lying within 200bp of TSS

gencode_v12_gtf <- file.path('~/projects/su2c_v2/data/raw/gtf_files/gencode.v12.annotation.gtf')
gencode_v12_txdb <- GenomicFeatures::makeTxDbFromGFF(gencode_v12_gtf, format = 'gtf')
gencode_v12_transcripts <- GenomicFeatures::transcripts(gencode_v12_txdb)
names(gencode_v12_transcripts) <- gencode_v12_transcripts$tx_name
stopifnot(!any(duplicated(names(gencode_v12_transcripts))))
promoter_prep <- SummExpDR::create_granges_prep(gencode_v12_transcripts, operation = 'upstream', analysis_name = 'tss200')

illumina_probe_anno <- as.data.frame(minfi::getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))

illumina_probe_anno_su2c <- SummExpDR::annotate_probes(illumina_anno = illumina_probe_anno,
                                                        granges_prep_list = list(promoter = promoter_prep),
                                                        prefix = 'su2c',
                                                        n_cores = 8)

gene_probe_map_tss200 <- SummExpDR::make_gene_probe_mapping(illumina_probe_anno_su2c,
                                                             mapping_col = 'su2c_tss200',
                                                             n_cores = 8)

illumina_EPIC_hg19_10b4_gene_2_probe_tss200 <- gene_probe_map_tss200
save(illumina_EPIC_hg19_10b4_gene_2_probe_tss200, file = './data/illumina_EPIC_hg19_10b4_gene_2_probe_tss200.RData')
sessionInfo()
