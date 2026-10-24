main <- function() {
  args <- commandArgs(trailingOnly = TRUE); stopifnot(length(args) == 1L)
  .libPaths(c('@R_LIB@', .libPaths()))
  suppressPackageStartupMessages({library(Seurat); library(Matrix); library(SingleCellExperiment);
    library(scDblFinder); library(celda); library(BiocParallel); library(jsonlite); library(ggplot2)})
  project <- '@PROJECT@'
  root <- file.path(project, '00.corrected_analysis/GSE332623/qc')
  id <- args[1]; out <- file.path(root, id); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  manifest <- read.delim(file.path(project, '00.corrected_analysis/GSE332623/capture_manifest.tsv'), check.names = FALSE)
  entry <- manifest[manifest$library_id == id, , drop = FALSE]; stopifnot(nrow(entry) == 1L)
  seed <- 332623L + match(id, sort(manifest$library_id))
  code <- normalizePath(sub('^--file=', '', grep('^--file=', commandArgs(FALSE), value = TRUE)[1]))
  code_sha <- digest::digest(code, algo = 'sha256', file = TRUE)
  state <- function(stage, extra = list()) {
    value <- c(list(status = 'running', stage = stage, library_id = id, updated = format(Sys.time(), '%FT%T%z'),
                   source_sha256 = entry$input_sha256, qc_code_sha256 = code_sha, seed = seed), extra)
    tmp <- tempfile(tmpdir = out); write_json(value, tmp, auto_unbox = TRUE, pretty = TRUE, null = 'null', digits = NA)
    stopifnot(file.rename(tmp, file.path(out, 'progress.json')))
    cat(format(Sys.time(), '%FT%T%z'), id, stage, '\n'); flush.console()
  }
  atomic_rds <- function(value, path) {tmp <- tempfile(tmpdir = out); saveRDS(value, tmp, compress = FALSE); stopifnot(file.rename(tmp, path))}
  tryCatch({
    state('loading_corrected_RNA_ADT_counts')
    stopifnot(digest::digest(entry$input_rds, algo = 'sha256', file = TRUE) == entry$input_sha256)
    s <- readRDS(entry$input_rds)
    r <- GetAssayData(s, assay = 'RNA', layer = 'counts'); adt <- GetAssayData(s, assay = 'ADT', layer = 'counts')
    stopifnot(identical(colnames(r), colnames(adt)), ncol(r) == entry$source_barcodes)
    fm <- s[['RNA']][[]]; symbols <- as.character(fm$source_symbol)
    stopifnot(identical(rownames(fm), rownames(r)))
    md <- s[[]]; stopifnot(identical(rownames(md), colnames(r)))
    rm(s); gc()
    metrics_path <- file.path(out, '01_metrics.rds')
    if (file.exists(metrics_path)) {
      checkpoint <- readRDS(metrics_path)
      stopifnot(checkpoint$source_sha256 == entry$input_sha256, checkpoint$policy_version == 'RNA-QC-v1')
      df <- checkpoint$metrics; threshold <- checkpoint$threshold
      stopifnot(identical(df$cell_id, colnames(r)))
    } else {
      umi <- Matrix::colSums(r); genes <- Matrix::colSums(r > 0)
      fraction <- function(index) if (any(index)) 100 * Matrix::colSums(r[index, , drop = FALSE]) / pmax(umi, 1) else rep(0, length(umi))
      df <- data.frame(cell_id = colnames(r), nCount_RNA = umi, nFeature_RNA = genes,
                       percent_mt = fraction(grepl('^MT-', symbols)),
                       percent_ribo = fraction(grepl('^RP[SL][0-9]', symbols)),
                       percent_hb = fraction(grepl('^HB[ABDEGMQZ]', symbols)),
                       nCount_ADT = Matrix::colSums(adt), nFeature_ADT = Matrix::colSums(adt > 0),
                       stringsAsFactors = FALSE, row.names = colnames(r))
      limits <- function(x, lower_floor) {
        z <- log10(x + 1); center <- median(z); spread <- mad(z)
        c(lower = max(lower_floor, floor(10^(center - 3 * spread) - 1)),
          upper = max(ceiling(quantile(x, .999)), ceiling(10^(center + 4 * spread) - 1)))
      }
      ul <- limits(umi, 200); gl <- limits(genes, 100)
      threshold <- list(UMI_min = unname(ul[1]), UMI_max = unname(ul[2]), genes_min = unname(gl[1]),
                        genes_max = unname(gl[2]), mt_max = max(15, min(30, median(df$percent_mt) + 4 * mad(df$percent_mt))),
                        rule = 'Lower log10 3MAD outliers with permissive UMI200/genes100 floors; upper max(log10 4MAD,99.9th percentile); MT max(15%,median+4MAD), capped30%.',
                        rationale = 'Preserve both low-RNA lymphocyte and higher-complexity modes; do not remove the lower mode at a histogram valley. Ribo/HB/ADT are review signals only.')
      df$qc_low_umi <- df$nCount_RNA < threshold$UMI_min; df$qc_high_umi <- df$nCount_RNA > threshold$UMI_max
      df$qc_low_genes <- df$nFeature_RNA < threshold$genes_min; df$qc_high_genes <- df$nFeature_RNA > threshold$genes_max
      df$qc_high_mt <- df$percent_mt > threshold$mt_max
      df$passes_basic_qc <- !(df$qc_low_umi | df$qc_high_umi | df$qc_low_genes | df$qc_high_genes | df$qc_high_mt)
      atomic_rds(list(source_sha256 = entry$input_sha256, policy_version = 'RNA-QC-v1', threshold = threshold, metrics = df), metrics_path)
      write_json(threshold, file.path(out, 'thresholds.json'), auto_unbox = TRUE, pretty = TRUE, digits = NA)
      for (metric in c('nCount_RNA','nFeature_RNA','percent_mt')) {
        vals <- df[[metric]]; d <- data.frame(value = vals)
        boundary <- if (metric == 'nCount_RNA') c(threshold$UMI_min, threshold$UMI_max) else
                    if (metric == 'nFeature_RNA') c(threshold$genes_min, threshold$genes_max) else threshold$mt_max
        if (metric != 'percent_mt') {d$value <- log10(d$value + 1); boundary <- log10(boundary + 1)}
        p <- ggplot(d, aes(x = value)) + geom_histogram(bins = 100, fill = '#297f86', color = 'white', linewidth = .1) +
          geom_vline(xintercept = boundary, color = '#b44932', linetype = 2) + theme_classic(base_size = 12) +
          labs(title = id, x = if (metric == 'percent_mt') 'Mitochondrial RNA (%)' else paste0('log10(', metric, ' + 1)'), y = 'Source barcodes')
        ggsave(file.path(out, paste0('QC_', metric, '_thresholds.png')), p, width = 7, height = 4, dpi = 150)
      }
    }
    state('basic_RNA_QC_complete', list(n_input = nrow(df), n_basic_retained = sum(df$passes_basic_qc)))
    basic <- which(df$passes_basic_qc); stopifnot(length(basic) >= 100L)
    doublet_path <- file.path(out, '02_doublets.rds')
    if (file.exists(doublet_path)) {
      dbl <- readRDS(doublet_path); stopifnot(dbl$source_sha256 == entry$input_sha256,
                                             identical(dbl$cell_ids, df$cell_id[basic]))
    } else {
      state('scDblFinder_running', list(n_cells = length(basic)))
      RNGkind("L'Ecuyer-CMRG"); set.seed(seed)
      register(SerialParam(RNGseed = seed))
      partitions <- if (length(basic) > 50000L) {
        assignment <- sample(rep(seq_len(ceiling(length(basic) / 25000L)), length.out = length(basic)))
        paste0(id, '__part', assignment)
      } else rep(id, length(basic))
      sce <- SingleCellExperiment(assays = list(counts = r[, basic, drop = FALSE]))
      sce <- scDblFinder(sce, samples = partitions, clusters = TRUE, nfeatures = 1352L, dims = 20L,
                        BPPARAM = SerialParam(RNGseed = seed), verbose = TRUE)
      data <- as.data.frame(colData(sce)); stopifnot(identical(rownames(data), df$cell_id[basic]))
      dbl <- list(source_sha256 = entry$input_sha256, seed = seed, cell_ids = df$cell_id[basic],
                  partitions = partitions, data = data,
                  parameters = list(method = 'scDblFinder', nfeatures = 1352L, dims = 20L,
                    expected_doublet_rate = 'Package default per QC library; no donor/group labels used',
                    large_library_split = '>50000 basic-pass barcodes split deterministically to <=25000 per partition',
                    physical_capture_verified = FALSE, partition_basis = 'one published MEX library'))
      atomic_rds(dbl, doublet_path); rm(sce); gc()
    }
    df$doublet_score <- NA_real_; df$doublet_class <- 'not_tested_basic_QC_failed'; df$doublet_partition <- NA_character_
    df$doublet_score[basic] <- dbl$data$scDblFinder.score
    df$doublet_class[basic] <- as.character(dbl$data$scDblFinder.class)
    df$doublet_partition[basic] <- dbl$partitions
    stopifnot(all(df$doublet_class[basic] %in% c('singlet','doublet')), all(is.finite(df$doublet_score[basic])))
    state('scDblFinder_complete', list(n_doublets = sum(df$doublet_class == 'doublet')))
    ambient_path <- file.path(out, '03_ambient.rds')
    if (file.exists(ambient_path)) {
      ambient <- readRDS(ambient_path); stopifnot(ambient$source_sha256 == entry$input_sha256,
                                                 identical(ambient$cell_ids, df$cell_id[basic]))
    } else {
      state('decontX_running', list(n_cells = length(basic)))
      set.seed(seed + 1L)
      nonzero <- which(Matrix::rowSums(r[, basic, drop = FALSE]) > 0)
      sce <- SingleCellExperiment(assays = list(counts = r[nonzero, basic, drop = FALSE]))
      sce <- celda::decontX(sce, maxIter = 500L, seed = seed + 1L, verbose = TRUE)
      stopifnot(identical(colnames(sce), df$cell_id[basic]))
      contamination <- as.numeric(colData(sce)$decontX_contamination)
      corrected <- assay(sce, 'decontXcounts')
      stopifnot(length(contamination) == length(basic), all(is.finite(contamination)),
                all(contamination >= 0 & contamination <= 1), identical(colnames(corrected), df$cell_id[basic]))
      ambient <- list(source_sha256 = entry$input_sha256, seed = seed + 1L, cell_ids = df$cell_id[basic],
                      contamination = contamination, corrected_counts = corrected,
                      cell_data = as.data.frame(colData(sce)),
                      parameters = list(method = 'decontX', maxIter = 500L, background = 'No empty-droplet matrix supplied; estimated from cell clusters',
                                        raw_counts_preserved = TRUE, corrected_matrix_is_fractional = TRUE,
                                        ambient_filter = 'Flag/review only; not an automatic exclusion'))
      atomic_rds(ambient, ambient_path); rm(sce, corrected); gc()
    }
    df$ambient_contamination <- NA_real_; df$ambient_contamination[basic] <- ambient$contamination
    acut <- max(.5, min(.9, median(ambient$contamination) + 4 * mad(ambient$contamination)))
    df$ambient_review <- FALSE; df$ambient_review[basic] <- ambient$contamination > acut
    df$passes_QC <- df$passes_basic_qc & df$doublet_class == 'singlet'
    df$library_id <- id; df$capture_id <- as.character(entry$capture_id); df$pool_id <- as.character(entry$pool_id)
    df$source_rna_gsm <- as.character(entry$rna_gsm); df$source_adt_gsm <- as.character(entry$adt_gsm)
    df$source_condition <- 'Unresolved'; df$Group <- 'Unresolved'
    df$source_barcode <- as.character(md$source_barcode)
    stopifnot(!anyNA(df$source_barcode), length(df$source_barcode) == ncol(r))
    df$excluded_reason <- ifelse(!df$passes_basic_qc, 'RNA_metric_outlier',
                                ifelse(df$doublet_class == 'doublet', 'scDblFinder_doublet', 'retained'))
    atomic_rds(df, file.path(out, '04_cell_QC.rds'))
    con <- gzfile(file.path(out, 'per_cell_QC.tsv.gz'), 'wt')
    write.table(df, con, sep = '\t', row.names = FALSE, quote = TRUE); close(con)
    stats <- list(status = 'QC_stages_complete_awaiting_bridge_validation', library_id = id,
                  source_sha256 = entry$input_sha256, qc_code_sha256 = code_sha, seed = seed,
                  n_input = nrow(df), n_basic_pass = sum(df$passes_basic_qc),
                  n_doublets = sum(df$doublet_class == 'doublet'), n_retained = sum(df$passes_QC),
                  n_ambient_review_retained = sum(df$passes_QC & df$ambient_review),
                  ambient_review_cutoff = acut, ambient_quantiles = quantile(ambient$contamination, c(0,.1,.25,.5,.75,.9,.95,.99,1)),
                  ambient_is_filter = FALSE, clinical_status = 'Unresolved', thresholds = threshold,
                  versions = lapply(c('Seurat','Matrix','scDblFinder','celda','SingleCellExperiment','BiocParallel'),
                                     function(p) list(package = p, version = as.character(packageVersion(p)))))
    write_json(stats, file.path(out, 'qc_stages_receipt.json'), auto_unbox = TRUE, pretty = TRUE, digits = NA)
    p <- ggplot(df[basic, ], aes(x = ambient_contamination, fill = doublet_class)) +
      geom_histogram(bins = 80, position = 'identity', alpha = .5) + geom_vline(xintercept = acut, linetype = 2) +
      theme_classic(base_size = 12) + labs(title = id, x = 'decontX estimated contamination (review only)', y = 'Barcodes')
    ggsave(file.path(out, 'QC_ambient_review.png'), p, width = 7, height = 4, dpi = 150)
    state('QC_stages_complete_awaiting_bridge_validation', stats)
  }, error = function(e) {
    write_json(list(status = 'failed', library_id = id, error = conditionMessage(e),
                    updated = format(Sys.time(), '%FT%T%z')), file.path(out, 'failure.json'), auto_unbox = TRUE, pretty = TRUE)
    stop(e)
  })
}
main()
