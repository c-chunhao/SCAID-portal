## GSE332623: read the two measured modalities separately, without QC or clinical inference.
## Source MEX matrices mix Gene Expression and Antibody Capture rows.
## Library/pool identifiers are technical labels, not donor or disease assignments.
ReadGSE332623Library <- function(input_dir, library_id, inventory = NULL) {
  stopifnot(length(input_dir) == 1L, length(library_id) == 1L)
  if (!grepl('^GSE332623_TotalseqB_NEA[1-7]_ADT[1-4]$', library_id))
    stop('Unexpected GSE332623 library identifier: ', library_id)
  prefix <- file.path(input_dir, library_id)
  paths <- c(matrix = paste0(prefix, '_matrix.mtx.gz'),
             features = paste0(prefix, '_features.tsv.gz'),
             barcodes = paste0(prefix, '_barcodes.tsv.gz'))
  if (!all(file.exists(paths))) stop('Missing MEX member for ', library_id)
  read_table <- function(path) {
    con <- gzfile(path, 'rt'); on.exit(close(con))
    utils::read.delim(con, header = FALSE, sep = '\t', quote = '',
                     comment.char = '', stringsAsFactors = FALSE, check.names = FALSE)
  }
  features <- read_table(paths[['features']])
  barcodes <- read_table(paths[['barcodes']])[[1L]]
  if (ncol(features) != 3L || anyNA(features) || anyNA(barcodes) ||
      anyDuplicated(barcodes) || any(!nzchar(barcodes)))
    stop('Invalid feature table or barcode axis: ', library_id)
  names(features) <- c('source_feature_id', 'source_symbol', 'feature_type')
  if (!setequal(unique(features$feature_type), c('Gene Expression', 'Antibody Capture')))
    stop('Expected the two explicit RNA/ADT feature types: ', library_id)
  features$source_row <- seq_len(nrow(features))
  con <- gzfile(paths[['matrix']], 'rt'); on.exit(close(con), add = TRUE)
  mat <- methods::as(Matrix::readMM(con), 'CsparseMatrix')
  if (!identical(as.integer(dim(mat)), as.integer(c(nrow(features), length(barcodes)))))
    stop('MEX axes do not match feature/barcode files: ', library_id)
  if (any(!is.finite(mat@x)) || any(mat@x < 0) || any(mat@x != round(mat@x)))
    stop('Source count matrix is not finite, nonnegative integer-valued: ', library_id)
  colnames(mat) <- paste0(library_id, '_', barcodes)
  create_assay <- function(feature_type, assay) {
    selected <- which(features$feature_type == feature_type)
    meta <- features[selected, , drop = FALSE]
    # Match Seurat's feature-name sanitization, and retain the unmodified names/IDs in metadata.
    symbols <- make.unique(gsub('[_|]', '-', meta$source_symbol))
    counts <- mat[selected, , drop = FALSE]; rownames(counts) <- symbols
    meta$assay_feature_name <- symbols; rownames(meta) <- symbols
    list(counts = counts, metadata = meta)
  }
  rna <- create_assay('Gene Expression', 'RNA')
  adt <- create_assay('Antibody Capture', 'ADT')
  obj <- SeuratObject::CreateSeuratObject(counts = rna$counts, project = 'GSE332623',
                                         assay = 'RNA', min.cells = 0, min.features = 0)
  obj[['ADT']] <- SeuratObject::CreateAssay5Object(counts = adt$counts, min.cells = 0, min.features = 0)
  obj[['RNA']] <- SeuratObject::AddMetaData(obj[['RNA']], metadata = rna$metadata)
  obj[['ADT']] <- SeuratObject::AddMetaData(obj[['ADT']], metadata = adt$metadata)
  source_gsm <- NA_character_
  if (!is.null(inventory)) {
    required <- c('gse', 'gsm', 'sample_title')
    if (!all(required %in% colnames(inventory))) stop('Inventory schema missing')
    hits <- unique(inventory$gsm[inventory$gse == 'GSE332623' &
                   inventory$sample_title == sub('^GSE332623_', '', library_id)])
    hits <- hits[!is.na(hits) & grepl('^GSM[0-9]+$', hits)]
    if (length(hits) != 1L) stop('No unique explicit library-to-GSM mapping: ', library_id)
    source_gsm <- hits[[1L]]
  }
  obj$orig.ident <- library_id
  obj$Sample <- library_id
  obj$source_library_id <- library_id
  obj$source_pool <- sub('^.*_(NEA[1-7])_.*$', '\\1', library_id)
  obj$source_barcode <- barcodes
  obj$source_gsm <- source_gsm
  obj$GSM <- source_gsm
  obj$Dataset <- 'GSE332623'
  obj$Tissue <- 'PBMC'
  obj$Group <- 'Unresolved'
  obj$source_condition <- 'Unresolved'
  obj$source_cohort <- 'AD and matched healthy controls; mixed pooled libraries'
  obj$Disease_detail <- 'Unresolved'
  obj$Disease_subtype <- 'Unresolved'
  obj$donor_id <- NA_character_
  obj$age_group <- NA_character_
  obj$condition_assignment_status <- 'No verified cell-to-donor/clinical mapping available'
  obj@misc$source_input <- list(library_id = library_id, paths = as.list(paths),
       source_feature_rows = nrow(features), source_barcode_count = length(barcodes),
       source_nnzero = length(mat@x), feature_type_counts = as.list(table(features$feature_type)),
       processing = 'Counts-only; no QC, normalization, embedding, annotation or disease assignment',
       source_study = 'GSE332623; DOI 10.1186/s12967-026-08288-7')
  stopifnot(identical(colnames(obj), paste0(library_id, '_', barcodes)),
            identical(colnames(obj[['RNA']]), colnames(obj[['ADT']])),
            nrow(obj[['RNA']]) == sum(features$feature_type == 'Gene Expression'),
            nrow(obj[['ADT']]) == sum(features$feature_type == 'Antibody Capture'),
            all(obj[['RNA']][[]]$feature_type == 'Gene Expression'),
            all(obj[['ADT']][[]]$feature_type == 'Antibody Capture'))
  obj
}

ReadGSE332623Libraries <- function(input_dir, inventory = NULL, library_ids = NULL) {
  if (is.null(library_ids)) {
    files <- list.files(input_dir, pattern = '^GSE332623_TotalseqB_NEA[1-7]_ADT[1-4]_matrix\\.mtx\\.gz$', full.names = FALSE)
    library_ids <- sort(sub('_matrix\\.mtx\\.gz$', '', files))
  }
  if (!length(library_ids) || anyDuplicated(library_ids)) stop('Invalid library list')
  objs <- setNames(vector('list', length(library_ids)), library_ids)
  for (u in library_ids) objs[[u]] <- ReadGSE332623Library(input_dir, u, inventory)
  attr(objs, 'skipped') <- character(0)
  objs
}
