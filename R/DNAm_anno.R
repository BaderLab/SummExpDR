

## Functions for generating DNA methylation annotations

#' Check if probes within within specified distance upstream of coordinate
#' @param coord = query coordinate
#' @param chr = chromosome indexed
#' @param regions = GenomicRanges object
#' @param bp = filter for probes within this distance of TSS
#' @return logical vector indicating for which genes the coord is upstream of
#' @export

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
#' @export
.check_isin <- function(coord, chr, regions) {
  range_vals <- GenomicRanges::ranges(regions)
  starts <- as.vector(GenomicRanges::start(range_vals), mode = 'integer')
  ends <- as.vector(GenomicRanges::end(range_vals), mode = 'integer')
  chr_annot <- as.vector(GenomicRanges::seqnames(regions), mode = 'character')
  return((coord >= starts & coord <= ends) & (chr == chr_annot))
}

#' check validity of granges_prep object
#' @param granges_prep granges_prep object
#' @export

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
#' @export

setClass('granges_prep',
         slots = c(granges_obj = 'GRanges',
                   operation = 'character',
                   analysis_name = 'character',
                   bp = 'integer'),
         validity = check_granges_prep)



#' create granges_prep
#' @description create 
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
#' @export
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
#' @return modified illumina_anno data frame with columns added for gene body annotations (semicolon delimited list of genes, '' if empty set),
#' promoter annotations (semicolon delimited list of genes, '' if empty set), and enhancer annotations (semicolon delimited list of start and end points
#' of enhancers). 
#' @export
annotate_probes <- function(illumina_anno, granges_prep_list, prefix, coord_column = 'pos', chr_column = 'chr', n_cores = 1) {
  stopifnot(coord_column %in% colnames(illumina_anno))
  stopifnot(chr_column %in% colnames(illumina_anno))
  stopifnot(is.integer(illumina_anno[,coord_column]))
  stopifnot(is.character(illumina_anno[,chr_column]))
  # run helper function in parallel or serial
  illumina_coords <- illumina_anno[,c(coord_column, chr_column)]
  n_probes <- nrow(illumina_coords)
  if (n_cores > 1) {
    # TODO: replace foreach implementation with BiocParallel implementation
    # # do parallel runs
    # cl <- snow::makeCluster(n_cores)
    # doSNOW::registerDoSNOW(cl)
    # pb <- txtProgressBar(max = n_probes, style = 3)
    # progress <- function(n) setTxtProgressBar(pb, n)
    # opts <- list(progress = progress)
    # `%dopar%` <- foreach::`%dopar%`
    # result_df <- foreach::foreach(i = 1:n_probes, .options.snow = opts, .export = c('.annotate_1_probe'), .combine = 'rbind') %dopar%
    # {
    #   return(.annotate_1_probe(illumina_coords = illumina_coords, 
    #                            probe_idx = i,
    #                            granges_prep_list = granges_prep_list))
    # }
    # snow::stopCluster(cl)
    
    
    
  } else {
    # do serial
    result_df <- c()
    pb <- txtProgressBar(max = nrow(illumina_coords), style = 3)
    for (i in 1:n_probes) {
      df_i <- .annotate_1_probe(illumina_coords = illumina_coords, 
                                probe_idx = i, 
                                granges_prep_list = granges_prep_list)
      result_df <- rbind(result_df, df_i)
      setTxtProgressBar(pb, value = i)
    }
    close(pb)
  }
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
#' @export

make_gene_probe_mapping <- function(illumina_anno,
                                    mapping_col,
                                    n_cores = 1) {
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
    # do parallel runs
    cl <- snow::makeCluster(n_cores)
    doSNOW::registerDoSNOW(cl)
    pb <- txtProgressBar(max = num_genes, style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    `%dopar%` <- foreach::`%dopar%`
    gene_probe_mapping_vect <- foreach::foreach(i = 1:num_genes, 
                                                .options.snow = opts, 
                                                .export = c('.map_1_gene_2_probes'), 
                                                .combine = 'c') %dopar% {
      gene_i <- genes_use[i]
      mapped_probes <- .map_1_gene_2_probes(mapping_vect, gene_i)
      return(mapped_probes)
    }
    snow::stopCluster(cl)
  } else {
    gene_probe_mapping_vect <- character(num_genes)
    pb <- txtProgressBar(min = 0, max = num_genes, style = 3)
    for (i in 1:num_genes) {
      gene_i <- genes_use[i]
      gene_probe_mapping_vect[i] <- .map_1_gene_2_probes(mapping_vect, gene_i)
      setTxtProgressBar(pb, i)
    }
    close(pb)
  }
  df_return <- data.frame(gene = genes_use, mapped_probes = gene_probe_mapping_vect, stringsAsFactors = FALSE)
  rownames(df_return) <- df_return$gene
  return(df_return)
}

.make_gmt_helper <- function(go_id,
                             biomart_anno,
                             gene_probe_mapping,
                             mapping_col,
                             biomart_gene_col = 'ensembl_transcript_id',
                             biomart_go_id_col = 'go_id',
                             biomart_go_type_col = 'namespace_1003',
                             biomart_go_term_col = 'name_1006',
                             min_genes = 15, 
                             max_genes = 200) {
  # go_id = GO id (character of length 1)
  # file_str = string to add 
  row_ind <- biomart_anno[,biomart_go_id_col] == go_id
  gene_ids_i <- unique(biomart_anno[,biomart_gene_col][row_ind])
  # skip if too many or too few genes associated with pathway
  if (length(gene_ids_i) < min_genes || length(gene_ids_i) > max_genes) {
   line_write <- '' 
  } else {
    probes_i <- character(0)
    for (ens_transcr_id_j in gene_ids_i) {
      # map genes to probes
      probes_i <- union(probes_i, strsplit(gene_probe_mapping[ens_transcr_id_j, mapping_col], split = ';')[[1]])
    }
    if (length(probes_i) > 0) {
      # add go term + probes to file
      go_type <- unique(biomart_anno[row_ind, biomart_go_type_col])
      go_term <- unique(biomart_anno[row_ind, biomart_go_term_col])
      line_write <- paste(paste(paste(toupper(go_term), 
                                      toupper(go_type), 
                                      go_id, sep = '%'),
                                go_term,
                                sep = '\t'), 
                          paste(probes_i, collapse = '\t'), 
                          sep = '\t')
    } else {
      line_write <- ''
    }
  }
  return(line_write)
}

#' Make GMT file mapping pathways to probes
#' @param  biomart_anno = biomart annotations containing ESNT ids + GO ids
#' @param gene_probe_mapping = dataframe with genes in rows, probes in a particular column, separated by semicolon delimited string
#' @param mapping_col = column in gene_probe_mapping
#' @param biomart_gene_col = gene column containing genes in biomart annotation
#' @param biomart_go_id_col = column in biomart_anno containing go ids
#' @param biomart_go_term_col = column in biomart_anno containing go terms
#' @param biomart_go_type_col = column in biomart_anno containing go types
#' @param min_genes = minimum # genes associated with a pathway. pathways with fewer will not be included in output file
#' @param max_genes = maximum # genes associated with a pathway. pathways with more will not be included in output file
#' @param n_cores = number of cores to use
#' @param output_file = file path to the output gmt file
#' @return  integer. 0 for completion without error, 1 for error
#' @export

DNAm_make_gmt <- function(biomart_anno,
                          gene_probe_mapping,
                          mapping_col,
                          biomart_gene_col = 'ensembl_transcript_id',
                          biomart_go_id_col = 'go_id',
                          biomart_go_type_col = 'namespace_1003',
                          biomart_go_term_col = 'name_1006',
                          min_genes = 15, 
                          max_genes = 200,
                          n_cores = 1,
                          output_file = file.path(output_dir, 'DNAm_probes.gmt')) {
  
  # check the gene to probe mapping
  mapping_regex <- '(([A-z0-9]|_)+;*([A-z0-9]|_)*)*'
  tryCatch({
    stopifnot(all(grepl(mapping_regex, gene_probe_mapping[,mapping_col])))
  }, error = function(e) {
    stop('invalid mapping column. check that mapping column has semicolon delimited strings')
  })
  tryCatch({
    bm_cols_check <- c(biomart_gene_col, biomart_go_id_col, biomart_go_term_col, biomart_go_type_col)
    bad_bm_cols <- setdiff(bm_cols_check, colnames(biomart_anno))
    stopifnot(length(bad_bm_cols) < 1)
  }, error = function(e) {
    stop(print(paste('Following columns not in biomart_anno:', paste(bad_bm_cols))))
  })
  tryCatch({
    stopifnot(mapping_col %in% colnames(gene_probe_mapping))
  }, error = function(e) {
    stop('mapping_col value is not a column of gene_probe_mapping')
  })
  # get rid of all genes without associated probe annotation
  genes_remove <- gene_probe_mapping[,mapping_col] == ''
  if (any(genes_remove)) {
    print(paste(length(genes_remove), 'of', nrow(gene_probe_mapping), 'genes in gene_probe_mapping',
                'have no associated probes, will not be used in mapping GO terms to probes'))
    gene_probe_mapping <- gene_probe_mapping[-which(genes_remove),]
  }
  # only use rows corresponding to genes in our gene to probe mapping
  common_genes <- intersect(rownames(gene_probe_mapping), biomart_anno[,biomart_gene_col])
  print(paste(length(common_genes), 
              'genes common to', 
              length(unique(biomart_anno[,biomart_gene_col])), 
              'biomart genes and', 
              nrow(gene_probe_mapping), 
              'gene_probe_mapping genes'))
  tryCatch({
    stopifnot(length(common_genes) > 0)
  }, error = function(e) {
    stop('number of common genes between biomart annotation and gene/probe mapping is < 1')
  })
  biomart_anno <- biomart_anno[which(biomart_anno[,biomart_gene_col] %in% common_genes), ]
  # get rid of all GO annotations that are IEA, RCA, ND
  biomart_anno <- biomart_anno[which(!biomart_anno$go_linkage_type %in% c('IEA', 'RCA', 'ND')), ]
  # get rid of all empty GO annotations
  empty_go <- biomart_anno[,biomart_go_id_col] == ''
  if (any(empty_go)) {
    biomart_anno <- biomart_anno[!empty_go,]
  }
  
  print(paste(length(unique(biomart_anno$ensembl_transcript_id)), ' of ', 
              nrow(gene_probe_mapping), ' gene ids from illumina annotation kept in GO term mappings'))
  
  # write gmt file
  if (file.exists(output_file)) {
    system(paste('rm', output_file))
  }
  print('Mapping pathways to EPIC probes')
  # loop through go ids and write file
  unique_go_ids <- unique(biomart_anno[,biomart_go_id_col])
  n_go_ids <- length(unique_go_ids)
  exit_status <- 1L
  try({
    if (n_cores > 1) {
      # do parallel runs
      cl <- snow::makeCluster(n_cores)
      doSNOW::registerDoSNOW(cl)
      pb <- txtProgressBar(max = n_go_ids, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n)
      opts <- list(progress = progress)
      `%dopar%` <- foreach::`%dopar%`
      
      gmt_lines <- foreach::foreach(i = 1:n_go_ids, .options.snow = opts, .export = c('.make_gmt_helper'), .combine = 'c') %dopar% {
        go_id_i <- unique_go_ids[i]
        line_i <- .make_gmt_helper(go_id = go_id_i,
                                   biomart_anno = biomart_anno,
                                   gene_probe_mapping = gene_probe_mapping,
                                   mapping_col = mapping_col,
                                   biomart_gene_col = biomart_gene_col,
                                   biomart_go_id_col = biomart_go_id_col,
                                   biomart_go_type_col = biomart_go_type_col,
                                   biomart_go_term_col = biomart_go_term_col,
                                   min_genes = min_genes, 
                                   max_genes = max_genes)
        return(line_i)
      }
      snow::stopCluster(cl)
    } else {
      gmt_lines <- c()
      pb <- txtProgressBar(min = 0, max = n_go_ids, style = 3)
      for (i in 1:n_go_ids) {
        go_id_i <- unique_go_ids[i]
        line_i <- .make_gmt_helper(go_id = go_id_i,
                                   biomart_anno = biomart_anno,
                                   gene_probe_mapping = gene_probe_mapping,
                                   mapping_col = mapping_col,
                                   biomart_gene_col = biomart_gene_col,
                                   biomart_go_id_col = biomart_go_id_col,
                                   biomart_go_type_col = biomart_go_type_col,
                                   biomart_go_term_col = biomart_go_term_col,
                                   min_genes = min_genes, 
                                   max_genes = max_genes)
        gmt_lines <- c(gmt_lines, line_i)
        setTxtProgressBar(pb, i)
      }
    }
    lines_remove <- which(gmt_lines == '')
    if (length(lines_remove) > 0) {
      print(paste('removing', length(lines_remove), 'pathways with no annotated probes'))
      gmt_lines <- gmt_lines[-lines_remove]
    }
    file_str <- paste(gmt_lines, collapse = '\n')
    print(paste('writing gmt file', output_file))
    output_fileconn <- file(output_file)
    writeLines(file_str, con = output_fileconn)
    print(paste('Finished writing gmt file with', length(gmt_lines), 'pathways'))
    close(output_fileconn)
    exit_status <- 0L
  })
  return(exit_status)
}
