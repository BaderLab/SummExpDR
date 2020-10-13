

## Functions for generating DNA methylation annotations

#' Check if probes within within specified distance upstream of coordinate
#' @param coord = query coordinate
#' @param chr = chromosome indexed
#' @param regions = GenomicRanges object
#' @param bp = filter for probes within this distance of TSS
#' @return logical vector indicating for which genes the coord is upstream of

.check_upstream <- function(coord, chr, regions, bp = 200) {
  range_vals <- GenomicRanges::ranges(regions)
  starts <- as.vector(GenomicRanges::start(range_vals), mode = 'integer')
  ends <- as.vector(GenomicRanges::end(range_vals), mode = 'integer')
  chr_annot <- as.vector(GenomicRanges::seqnames(regions), mode = 'character')
  strand <- as.vector(GenomicRanges::strand(regions), mode = 'character')
  neg_strand <- strand == '-'
  if (!all(strand[!neg_strand] == '+')) {
    stop('expect strand to be + or -')
  }
  result <- logical(length(range_vals))
  start_dist <- (starts - coord)
  end_dist <- (coord - ends)
  in_chr <- (chr_annot == chr)
  if (any(!neg_strand)) {
    # for positive strand, difference between forward start and coordinate
    result[!neg_strand] <- (((start_dist <= bp) & (start_dist > 0)) & in_chr)[!neg_strand]
  }
  if (any(neg_strand)) {
    # for negative strand, difference between forward end and coordinate
    result[neg_strand] <- (((end_dist <= bp) & (end_dist > 0)) & in_chr)[neg_strand]
  }
  return(result)
}

#' Check if a probe is within a set of genomic regions
#' @param coord = query coordinate
#' @param chr = chromosome indexed
#' @param regions = GRanges object with start + end points
#' @return logical vector indicating which regions the coordinate lies within
.check_isin <- function(coord, chr, regions) {
  range_vals <- GenomicRanges::ranges(regions)
  starts <- as.vector(GenomicRanges::start(range_vals), mode = 'integer')
  ends <- as.vector(GenomicRanges::end(range_vals), mode = 'integer')
  chr_annot <- as.vector(GenomicRanges::seqnames(regions), mode = 'character')
  return((coord >= starts & coord <= ends) & (chr == chr_annot))
}

#' check validity of granges_prep object
#' @param granges_prep granges_prep object

check_granges_prep <- function(granges_prep) {
  stopifnot(is(granges_prep@granges_obj, 'GRanges'))
  stopifnot(is.character(granges_prep@operation))
  stopifnot(length(granges_prep@operation) == 1L)
  stopifnot(granges_prep@operation %in% c('in', 'upstream'))
  stopifnot(is.character(granges_prep@analysis_name))
  stopifnot(length(granges_prep@analysis_name) == 1L)
  stopifnot(is.integer(granges_prep@bp))
  stopifnot(length(granges_prep@bp) == 1L)
  stopifnot(!is.null(names(granges_prep@granges_obj)))
  stopifnot(!any(is.na(names(granges_prep@granges_obj))))
  stopifnot(is.character(names(granges_prep@granges_obj)))
}

#' Container For GRanges Object + Operation
#' @slot granges_obj = GRanges object, with names set
#' @slot operation = operation to perform, 1 of in/upstream
#' @slot analysis_name = name for granges analysis
#' @slot bp = integer: number of base pairs upstream if checking upstream

setClass('granges_prep',
         slots = c(granges_obj = 'GRanges',
                   operation = 'character',
                   analysis_name = 'character',
                   bp = 'integer'),
         validity = check_granges_prep)



#' create granges_prep
#' @description create granges_prep object
#' @param granges_obj = GRanges object
#' @param operation = operation to perform (string, one of 'in' or 'upstream')
#' @param analysis_name = name for granges analysis
#' @param bp = number of basepairs upstream
#' @export

create_granges_prep <- function(granges_obj, operation, analysis_name, bp = 200L) {
  new_obj <- new('granges_prep',
                 granges_obj = granges_obj,
                 operation = operation,
                 analysis_name = analysis_name,
                 bp = bp)
  return(new_obj)
}

#' Helper function used for annotate_probes
#' @param illumina_coords = data frame containing probe coordinates (column 1) and corresponding chromosomes (column 2)
#' @param probe_idx = probe index (integer to index illumina_coords by row)
#' @param transcript_grange = GenomicRanges object with transcript coordinates
#' @param enh_grange = GenomicRanges objec with enhancer coordinates
#' @return 1 row dataframe with a column for each granges_prep object in granges_prep_list
.annotate_1_probe <- function(illumina_coords, probe_idx, granges_prep_list) {
  probe_coord <- illumina_coords[probe_idx, 1]
  probe_chr <- illumina_coords[probe_idx, 2]
  out_list <- c()

  .annotate_1_granges <- function(probe_coord, probe_chr, granges_prep) {
    check_granges_prep(granges_prep)
    granges_obj <- granges_prep@granges_obj
    operation <- granges_prep@operation
    range_names <- names(granges_obj)
    if (operation == 'in') {
      is_within_range <- .check_isin(probe_coord, probe_chr, granges_obj)
    } else {
      bp <- granges_prep@bp
      is_within_range <- .check_upstream(probe_coord, probe_chr, granges_obj, bp = bp)
    }
    if (any(is_within_range)) {
      names_use <- range_names[which(is_within_range)]
      anno_str <- paste(names_use, collapse = ';')
    } else {
      anno_str <- ''
    }
    return(anno_str)
  }

  for (granges_prep in granges_prep_list) {
    anno_str <- .annotate_1_granges(probe_coord, probe_chr, granges_prep)
    new_list <- list(anno_str)
    names(new_list)[1] <- granges_prep@analysis_name
    out_list <- c(out_list, new_list)
    rm(new_list)
  }
  out_df <- data.frame(out_list, stringsAsFactors = FALSE)
  return(out_df)
}

#' build annotations
#' @param illumina_anno = data frame from illumina EPIC array manifest file, must have MAPINFO column containing cpg coordinates
#' @param granges_prep_list = list of granges_prep objects
#' @param coord_column = column for coordinates of methylation sites
#' @param chr_column = column with chromosomes that probes are annotated to
#' @param prefix = column name prefix for annotations
#' @param pbar = whether or not to add progress bar
#' @return modified illumina_anno data frame with columns added for gene body annotations (semicolon delimited list of genes, '' if empty set),
#' promoter annotations (semicolon delimited list of genes, '' if empty set), and enhancer annotations (semicolon delimited list of start and end points
#' of enhancers).
#' @export
annotate_probes <- function(illumina_anno, granges_prep_list, prefix, coord_column = 'pos', chr_column = 'chr', n_cores = 1, pbar = TRUE) {
  stopifnot(coord_column %in% colnames(illumina_anno))
  stopifnot(chr_column %in% colnames(illumina_anno))
  stopifnot(is.integer(illumina_anno[,coord_column]))
  stopifnot(is.character(illumina_anno[,chr_column]))
  # run helper function in parallel or serial
  probe_coords <- illumina_anno[,c(coord_column, chr_column)]
  n_probes <- nrow(probe_coords)
  if (n_cores > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  }
  df_list <- BiocParallel::bplapply(1:n_probes, BPPARAM = bpparam, FUN = function(x) {
    .annotate_1_probe(probe_idx = x, illumina_coords = probe_coords, granges_prep_list = granges_prep_list)
    })
  result_df <- as.data.frame(dplyr::bind_rows(df_list))
  # set column names for added annotations
  colnames(result_df) <- paste(prefix, colnames(result_df), sep = '_')
  illumina_anno <- cbind(illumina_anno, result_df)
  return(illumina_anno)
}

.map_1_gene_2_probes <- function(mapping_vect, gene) {
  # mapping_vect = mapping_vect in make_gene_probe_mapping
  # gene = character vector denoting gene of interest
  # returns: string denoting all probes that map to the gene in the mapping_vect (probe to gene mapping)
  gene_regex <- paste0('(^|;)', gene, '($|;)')
  matched_inds <- grep(gene_regex, mapping_vect)
  if (length(matched_inds) > 0) {
    probe_names <- names(mapping_vect)[matched_inds]
    probe_names <- paste(probe_names, collapse = ';')
  } else {
    # note that this case should not happen
    probe_names <- ''
  }
  return(probe_names)
}
#' Make gene to probe mapping
#' @param illumina_anno = data.frame or DataFrame of probe metadata, containing
#' a mapping of probes to genes for a particular type of annotation, e.g. gene body, enhancer, promoter
#' @param mapping_col = name of column containing mapping. column contains character vector of semicolon delimited
#' strings denoting genes each probe maps to
#' @param n_cores = number of cores to use
#' @param pbar = whether to include progress bar
#' @export

make_gene_probe_mapping <- function(illumina_anno,
                                    mapping_col,
                                    n_cores = 1,
                                    pbar = TRUE) {
  if (!mapping_col %in% colnames(illumina_anno)) {
    stop(paste(mapping_col, 'is not a column in', illumina_anno))
  }
  # get mapping vect and run check
  mapping_vect <- illumina_anno[,mapping_col]
  names(mapping_vect) <- rownames(illumina_anno)
  mapping_regex <- '(([A-z0-9]|_)+;*([A-z0-9]|_)*)*'
  tryCatch({
    stopifnot(all(grepl(mapping_regex, mapping_vect)))
  }, error = function(e) {
    stop('invalid mapping column. check that mapping column has semicolon delimited strings')
  })
  # remove empty indices
  empty_inds <- mapping_vect == ''
  if (any(empty_inds)) {
    mapping_vect <- mapping_vect[-which(empty_inds)]
  }
  # get all unique genes
  genes_use <- unique(unlist(strsplit(mapping_vect, split = ';')))
  if (any(genes_use == '')) {
    genes_use <- genes_use[-which(genes_use == '')]
  }
  num_genes <- length(genes_use)
  if (n_cores > 1) {
    bpparam <- BiocParallel::MulticoreParam(workers = n_cores, progressbar = pbar)
  } else {
    bpparam <- BiocParallel::SerialParam(progressbar = pbar)
  }
  mapping_list <- BiocParallel::bplapply(genes_use, BPPARAM = bpparam, FUN = function(x) {.map_1_gene_2_probes(mapping_vect, gene = x)})
  df_return <- data.frame(gene = genes_use, mapped_probes = unlist(mapping_list), stringsAsFactors = FALSE)
  rownames(df_return) <- df_return$gene
  return(df_return)
}
