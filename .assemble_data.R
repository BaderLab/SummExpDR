set.seed(42L)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(SummExpDR)
library(GenomicFeatures)
library(GEOquery)
library(minfi)
library(SummarizedExperiment)
library(biomaRt)
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

# assemble gmt file for extra test data dealing with set overlap
gmt_file <- '../su2c_v2/data/raw/gmt_files/baderlab_gmts/Human_GO_AllPathways_no_GO_iea_April_01_2018_symbol.gmt'
baderlab_human_pathways_dec_2018 <- SummExpDR::load_gmt(gmt_file)
save(baderlab_human_pathways_dec_2018, file = './data/baderlab_human_pathways_dec_2018.RData')

# Assemble RNA-seq data from GSE77938
tmp <- tempdir()
dir.create(tmp)
data_link <- 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77938/suppl/GSE77938_discovery_gene_counts.txt.gz'
base_name_counts <- basename(data_link)
curl::curl_download(data_link, destfile = file.path(tmp, base_name_counts))
system(paste('gunzip', file.path(tmp, base_name_counts)))
txt_file_counts <- sub('.gz$', '', base_name_counts)
counts_df <- read.delim(file.path(tmp, txt_file_counts), row.names = 1)
counts_mat <- as.matrix(counts_df)
counts_mat <- counts_mat[order(matrixStats::rowVars(counts_mat), decreasing = TRUE)[1:5000], ]
RNA_counts <- counts_mat

# make some made up metadata
RNA_meta <- S4Vectors::DataFrame(KC = grepl('KC', colnames(RNA_counts)),
                                 Age = c(rnorm(12, 54, sd = 20), rnorm(4, 75, 20)),
                                 Sex = c(rep('M', 5), rep('F', 8), rep('M', '3')),
                                 row.names = colnames(RNA_counts))

RNA_SE <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = RNA_counts), colData = RNA_meta)

save(RNA_SE, file =  'data/RNA_SE.RData')
system(paste('rm -rf', tmp))

# Assemble Microarray Data from Cavalli et al. 2017
tmp <- tempdir()
dir.create(tmp)
GSEM <- GEOquery::getGEO(GEO = 'GSE85217', GSEMatrix = T, destdir = tmp)
GSEM <- GSEM$GSE85217_series_matrix.txt.gz
metadata <- GSEM@phenoData@data
accession <- metadata$geo_accession
rownames(metadata) <- accession
print('Loading ')
for (i in 1:length(accession)) {
  id.i <- accession[i]
  sample.i <- GEOquery::getGEO(GEO = accession[i], GSEMatrix = T)
  stopifnot(is(sample.i, 'GSM'))
  genes <- sample.i@dataTable@table$ID_REF
  value <- sample.i@dataTable@table$VALUE
  if (i == 1) {
    num.genes <- length(genes)
    num.samples <- length(accession)
    new.mat <- matrix(numeric(num.genes*num.samples),
                      nrow = num.genes,
                      ncol = num.samples)
    rownames(new.mat) <- genes
    colnames(new.mat) <- accession
  }
  new.mat[genes ,id.i] <- value
}
cavalli_2017_mat <- new.mat

archives <- biomaRt::listEnsemblArchives()
# archives[grep('Ensembl 77', archives[,1]),]
# version                                 date
# "Ensembl 77"                           "Oct 2014"
# url
# "http://Oct2014.archive.ensembl.org"
ensembl <- biomaRt::useMart(biomart = 'ENSEMBL_MART_ENSEMBL', ensemblRedirect = F,
                   host = 'http://Oct2014.archive.ensembl.org')
# datasets <- listDatasets(ensembl)
# datasets[grep('hsapiens', datasets$dataset),]
# dataset                 description version
# 32 hsapiens_gene_ensembl Homo sapiens genes (GRCh38)  GRCh38
ensembl <- biomaRt::useDataset(mart = ensembl, dataset = 'hsapiens_gene_ensembl')
attr <- biomaRt::listAttributes(ensembl)
# attr[grep('hgnc', attr$name),]
# name          description         page
# 61              hgnc_id           HGNC ID(s) feature_page
# 62          hgnc_symbol          HGNC symbol feature_page
# 63 hgnc_transcript_name HGNC transcript name feature_page
# head(attr)
# name           description         page
# 1       ensembl_gene_id       Ensembl Gene ID feature_page
# 2 ensembl_transcript_id Ensembl Transcript ID feature_page
# 3    ensembl_peptide_id    Ensembl Protein ID feature_page
# 4       ensembl_exon_id       Ensembl Exon ID feature_page
# 5           description           Description feature_page
# 6       chromosome_name       Chromosome Name feature_page

cavalli.genes <- rownames(cavalli_2017_mat)
cavalli.genes <- sub('_at$', '', cavalli.genes)
rownames(cavalli_2017_mat) <- cavalli.genes
# head(cavalli.genes)
# [1] "ENSG00000000003" "ENSG00000000005" "ENSG00000000419" "ENSG00000000457" "ENSG00000000460"
# [6] "ENSG00000000938"
gene.mapping <- biomaRt::getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                              mart = ensembl,
                              filters = c('ensembl_gene_id'),
                              values = cavalli.genes)

# all(cavalli.genes %in% gene.mapping$ensembl_gene_id)
# [1] TRUE
# all(gene.mapping$ensembl_gene_id %in% cavalli.genes)
# [1] TRUE
# length(cavalli.genes)
# [1] 21641
gene.mapping <- gene.mapping[which(!duplicated(gene.mapping$hgnc_symbol)),]
gene.mapping <- gene.mapping[which(!duplicated(gene.mapping$ensembl_gene_id)),]
# '' %in% gene.mapping$hgnc_symbol
# # [1] TRUE
# gene.mapping[which(gene.mapping$hgnc_symbol == ''),]
# ensembl_gene_id hgnc_symbol
# 85 ENSG00000005189
gene.mapping <- gene.mapping[-which(gene.mapping$hgnc_symbol == ''),]
if (!all(grepl('[A-Z0-9]+', gene.mapping$hgnc_symbol))) {
  stop('expect all genes to match HGNC regex')
}
# [1] TRUE
cavalli.genes <- cavalli.genes[cavalli.genes %in% gene.mapping$ensembl_gene_id]
cavalli.genes <- cavalli.genes[order(cavalli.genes)]
# length(cavalli.genes)
# [1] 20223
rownames(gene.mapping) <- gene.mapping$ensembl_gene_id


# reorder/subset matrix, add hgnc symbols
cavalli.hgnc <- gene.mapping[cavalli.genes, 'hgnc_symbol']
cavalli_2017_mat <- cavalli_2017_mat[cavalli.genes,]
rownames(cavalli_2017_mat) <- cavalli.hgnc

if (length(intersect(rownames(metadata), colnames(cavalli_2017_mat))) != nrow(metadata)) {
  stop('expect identical samples in metadata and expression matrix')
} else if (!all(rownames(metadata) == colnames(cavalli_2017_mat))) {
  metadata <- metadata[colnames(cavalli_2017_mat),]
}

cavalli_SE_hgnc <- SummarizedExperiment::SummarizedExperiment(assays = list(RNA = cavalli_2017_mat),
                                                              colData = S4Vectors::DataFrame(metadata))

set.seed(12345)
all.samples <- colnames(cavalli_SE_hgnc)
subs.samples <- sample(x = all.samples, size = floor(0.25*length(all.samples)))
cavalli_2017_testdata <- cavalli_SE_hgnc[,subs.samples]
save(cavalli_2017_testdata, file = file.path('data', 'cavalli_2017_testdata.RData'))

system(paste('rm -rf', tmp))

sessionInfo()
