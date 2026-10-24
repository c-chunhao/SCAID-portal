main <- function() {
  args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1L)
  .libPaths(c('@R_LIB@', .libPaths()))
  suppressPackageStartupMessages({library(Seurat); library(Matrix); library(rhdf5); library(jsonlite)})
  project <- '@PROJECT@'; id <- args[1]
  out <- file.path(project, '00.corrected_analysis/GSE332623/qc', id)
  entry <- read.delim(file.path(project, '00.corrected_analysis/GSE332623/capture_manifest.tsv'), check.names = FALSE)
  entry <- entry[entry$library_id == id, , drop = FALSE]; stopifnot(nrow(entry) == 1L)
  stats <- read_json(file.path(out, 'qc_stages_receipt.json'), simplifyVector = TRUE)
  stopifnot(stats$status == 'QC_stages_complete_awaiting_bridge_validation',
            stats$source_sha256 == entry$input_sha256)
  target <- file.path(out, 'counts_bridge.h5'); rds_target <- file.path(out, 'filtered_rna_adt_counts.rds')
  if (file.exists(file.path(out, 'receipt.json'))) {
    previous <- read_json(file.path(out, 'receipt.json'), simplifyVector = TRUE)
    stopifnot(previous$status == 'validated', previous$source_sha256 == entry$input_sha256,
              digest::digest(target, algo = 'sha256', file = TRUE) == previous$bridge_sha256,
              digest::digest(rds_target, algo = 'sha256', file = TRUE) == previous$filtered_rds_sha256)
    cat('Validated QC output already present:', id, '\n'); return(invisible(NULL))
  }
  df <- readRDS(file.path(out, '04_cell_QC.rds')); ambient <- readRDS(file.path(out, '03_ambient.rds'))
  stopifnot(digest::digest(entry$input_rds, algo = 'sha256', file = TRUE) == entry$input_sha256)
  s <- readRDS(entry$input_rds); r <- GetAssayData(s, assay = 'RNA', layer = 'counts'); adt <- GetAssayData(s, assay = 'ADT', layer = 'counts')
  rf <- s[['RNA']][[]]; af <- s[['ADT']][[]]; original_md <- s[[]]
  stopifnot(identical(df$cell_id, colnames(r)), identical(colnames(r), colnames(adt)),
            identical(rownames(rf), rownames(r)), identical(rownames(af), rownames(adt)))
  keep <- which(df$passes_QC); stopifnot(length(keep) == stats$n_retained, length(keep) >= 100)
  cells <- df$cell_id[keep]
  RNA <- r[, keep, drop = FALSE]; ADT <- adt[, keep, drop = FALSE]
  ix <- match(cells, ambient$cell_ids); stopifnot(!anyNA(ix))
  DC <- ambient$corrected_counts[, ix, drop = FALSE]
  ## Pad genes dropped only for zero count during decontX fitting, without touching raw RNA.
  map <- match(rownames(DC), rownames(RNA)); stopifnot(!anyNA(map), !anyDuplicated(map), all(diff(map) > 0))
  DC@i <- as.integer(map[DC@i + 1L] - 1L); DC@Dim[1] <- nrow(RNA); DC@Dimnames <- dimnames(RNA)
  validObject(DC)
  md <- original_md[keep, , drop = FALSE]
  for (name in setdiff(colnames(df), 'cell_id')) md[[name]] <- df[[name]][keep]
  md$source_mex_gsm <- as.character(original_md$source_gsm[keep])
  md$source_gsm <- as.character(entry$rna_gsm)
  md$source_rna_gsm <- as.character(entry$rna_gsm); md$source_adt_gsm <- as.character(entry$adt_gsm)
  md$library_id <- id; md$capture_id <- as.character(entry$capture_id); md$pool_id <- as.character(entry$pool_id)
  md$physical_capture_verified <- FALSE; md$clinical_assignment <- 'Unresolved'
  stopifnot(identical(rownames(md), cells), all(md$passes_QC), all(md$Group == 'Unresolved'))
  rm(s, r, adt, original_md, ambient); gc()
  temp <- tempfile(pattern = '.bridge-', tmpdir = out, fileext = '.h5')
  on.exit(unlink(temp), add = TRUE); h5createFile(temp)
  fid <- H5Fopen(temp)
  h5writeAttribute('scaid-rna-adt-counts-v1', fid, 'schema', asScalar = TRUE,
                   encoding = 'UTF-8', variableLengthString = TRUE)
  H5Fclose(fid)
  vector_write <- function(path, value, type = NULL) {
    if (is.character(value)) {
      value[is.na(value)] <- ''; h5write(enc2utf8(value), temp, path, encoding = 'UTF-8', level = 1)
    } else {
      value <- if (is.logical(value)) as.integer(value) else value
      h5createDataset(temp, path, dims = length(value), chunk = max(1L, min(length(value), 1000000L)),
                      storage.mode = if (is.integer(value)) 'integer' else 'double', H5type = type, level = 1)
      h5write(value, temp, path)
    }
  }
  verify_vector <- function(path, expected) {
    for (start in seq.int(1L, length(expected), by = 1000000L)) {
      count <- min(1000000L, length(expected) - start + 1L)
      got <- h5read(temp, path, start = start, count = count)
      stopifnot(isTRUE(all.equal(as.vector(got), as.vector(expected[start:(start + count - 1L)]),
                                check.attributes = FALSE, tolerance = 0)))
    }
  }
  write_matrix <- function(name, m, features, raw = TRUE) {
    stopifnot(identical(colnames(m), cells), identical(rownames(features), rownames(m)),
              all(is.finite(m@x)), all(m@x >= 0))
    if (raw) stopifnot(all(m@x == floor(m@x)))
    h5createGroup(temp, name)
    vector_write(paste0(name, '/data'), m@x, 'H5T_IEEE_F64LE')
    vector_write(paste0(name, '/indices'), m@i, 'H5T_STD_I32LE')
    vector_write(paste0(name, '/indptr'), m@p, 'H5T_STD_I64LE')
    vector_write(paste0(name, '/shape'), as.integer(c(ncol(m), nrow(m))), 'H5T_STD_I64LE')
    vector_write(paste0(name, '/feature_ids'), as.character(features$source_feature_id))
    vector_write(paste0(name, '/feature_names'), rownames(m))
    vector_write(paste0(name, '/source_symbols'), as.character(features$source_symbol))
    vector_write(paste0(name, '/source_row'), as.integer(features$source_row))
    for (field in c('data','indices','indptr')) verify_vector(paste0(name, '/', field), switch(field, data = m@x, indices = m@i, indptr = m@p))
    stopifnot(identical(as.character(h5read(temp, paste0(name, '/feature_names'))), rownames(m)))
    list(shape = c(ncol(m), nrow(m)), nnz = length(m@x), sum = sum(m@x), source_integer_counts = raw)
  }
  matrix_records <- list(RNA = write_matrix('RNA', RNA, rf, TRUE), ADT = write_matrix('ADT', ADT, af, TRUE),
                         RNA_decontX = write_matrix('RNA_decontX', DC, rf, FALSE))
  h5createGroup(temp, 'obs'); vector_write('obs/cell_id', cells)
  for (name in colnames(md)) {
    stopifnot(!grepl('/', name, fixed = TRUE))
    value <- md[[name]]
    if (is.factor(value)) value <- as.character(value)
    vector_write(paste0('obs/', name), value)
  }
  stopifnot(identical(as.character(h5read(temp, 'obs/cell_id')), cells),
            identical(as.numeric(h5read(temp, 'RNA/shape')), as.numeric(c(length(cells), nrow(RNA)))),
            identical(as.numeric(h5read(temp, 'ADT/shape')), as.numeric(c(length(cells), nrow(ADT)))))
  rds_temp <- tempfile(pattern = '.filtered-', tmpdir = out, fileext = '.rds')
  on.exit(unlink(rds_temp), add = TRUE)
  saveRDS(list(RNA = RNA, ADT = ADT, RNA_decontX = DC, cell_metadata = md,
               rna_feature_metadata = rf, adt_feature_metadata = af,
               provenance = list(source_rds = entry$input_rds, source_sha256 = entry$input_sha256,
                 raw_counts_preserved = TRUE, RNA_decontX_is_fractional_processed_counts = TRUE,
                 clinical_status = 'Unresolved', qc_stages = stats)), rds_temp, compress = 'gzip')
  check <- readRDS(rds_temp)
  stopifnot(identical(check$RNA, RNA), identical(check$ADT, ADT), identical(check$RNA_decontX, DC),
            identical(check$cell_metadata, md)); rm(check); gc()
  if (file.exists(target) || file.exists(rds_target)) stop('Unvalidated existing output needs inspection; refusing overwrite')
  stopifnot(file.rename(temp, target), file.rename(rds_temp, rds_target))
  receipt <- c(list(status = 'validated', library_id = id, capture_id = as.character(entry$capture_id),
                   physical_capture_verified = FALSE, source_sha256 = entry$input_sha256,
                   bridge_path = target, bridge_sha256 = digest::digest(target, algo = 'sha256', file = TRUE),
                   filtered_rds_path = rds_target, filtered_rds_sha256 = digest::digest(rds_target, algo = 'sha256', file = TRUE),
                   cell_QC_table = file.path(out, 'per_cell_QC.tsv.gz'),
                   clinical_status = 'Unresolved', raw_RNA_ADT_preserved = TRUE,
                   validation = list(all_sparse_components_readback_equal = TRUE, raw_counts_integer = TRUE,
                     RNA_ADT_same_cell_axis = TRUE, metadata_axis_equal = TRUE, filtered_RDS_roundtrip = TRUE),
                   matrices = matrix_records, completed = format(Sys.time(), '%FT%T%z')),
               stats[setdiff(names(stats), c('status','library_id','source_sha256','clinical_status'))])
  tmp <- tempfile(tmpdir = out); write_json(receipt, tmp, auto_unbox = TRUE, pretty = TRUE, null = 'null', digits = NA)
  stopifnot(file.rename(tmp, file.path(out, 'receipt.json')))
  cat('QC bridge validated:', id, length(cells), 'retained cells\n')
}
main()
