# Scientific reproducibility/axis safeguards for versioned SCAID integration.
# These functions mutate only caller-owned objects. Historical release RDS files
# must be preserved; callers write a distinct vNext output path.

freeze_analysis_rng <- function(seed = 20261005L) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  options(mc.cores = 4L)
  invisible(list(seed=seed, RNGkind=RNGkind(), threads=4L))
}

source_cell_key <- function(meta) {
  need <- c('Dataset','Sample','OldCells')
  stopifnot(all(need %in% names(meta)))
  if(anyNA(meta[,need]) || any(meta[,need]=='') || any(vapply(meta[,need], function(z)any(grepl('|',z,fixed=TRUE)),logical(1))))
    stop('Incomplete or ambiguous source cell key')
  key <- do.call(paste,c(meta[,need],sep='|'))
  if(anyDuplicated(key)) stop('Duplicate exact source cells; reconcile before integration')
  key
}

get_RNA_counts_checked <- function(object) {
  if(!'RNA'%in%names(object@assays))stop('RNA assay missing; do not use integrated data as counts')
  if(inherits(object[['RNA']],'Assay5')) {
    layers<-SeuratObject::Layers(object[['RNA']],search='^counts')
    if(length(layers)!=1L||layers!='counts')stop('RNA has split counts layers; explicitly join/iterate them with cell provenance before retrieving counts')
  }
  counts <- if(utils::packageVersion('SeuratObject') >= '5.0.0')
    SeuratObject::GetAssayData(object,assay='RNA',layer='counts') else
    SeuratObject::GetAssayData(object,assay='RNA',slot='counts')
  if(!inherits(counts,'sparseMatrix'))stop('Expected sparse unintegrated RNA counts')
  if(!identical(colnames(counts),SeuratObject::Cells(object)))stop('RNA counts cell axis differs from object')
  counts
}

attach_reduction_checked <- function(object, embeddings, name, assay='RNA') {
  cells <- SeuratObject::Cells(object)
  if(!is.matrix(embeddings) || is.null(rownames(embeddings)) || anyDuplicated(rownames(embeddings)))
    stop('Embedding requires unique, explicit cell IDs')
  idx <- match(cells,rownames(embeddings))
  if(anyNA(idx)) stop('Embedding misses cells; do not attach positional rows')
  aligned <- embeddings[idx,,drop=FALSE]
  if(!identical(rownames(aligned),cells) || !all(is.finite(aligned))) stop('Invalid aligned embedding')
  object[[name]] <- SeuratObject::CreateDimReducObject(embeddings=aligned, key=paste0(gsub('[^A-Za-z0-9]','',name),'_'),assay=assay)
  object
}

assemble_annotations_checked <- function(base_object, metadata_list,
                                         columns=c('celltype1','celltype2')) {
  stopifnot(is.list(metadata_list),length(metadata_list)>0)
  metadata_list <- lapply(metadata_list,function(m){
    if(!all(columns%in%names(m)) || is.null(rownames(m))) stop('Each annotation layer needs named rows and explicit columns')
    m[,columns,drop=FALSE]
  })
  md <- do.call(rbind,metadata_list)
  if(anyDuplicated(rownames(md))) stop('Annotation branches contain duplicate cell IDs')
  if(anyNA(md)||any(as.matrix(md)=='')) stop('Incomplete annotation layers')
  if(!all(rownames(md)%in%SeuratObject::Cells(base_object))) stop('Annotation cells absent from the explicitly named base object')
  x <- subset(base_object,cells=rownames(md))
  md <- md[match(SeuratObject::Cells(x),rownames(md)),,drop=FALSE]
  stopifnot(identical(rownames(md),SeuratObject::Cells(x)))
  SeuratObject::AddMetaData(x,metadata=md)
}

full_RNA_log_normalization <- function(counts,target_sum=10000) {
  if(!inherits(counts,'sparseMatrix')) stop('Use sparse full-RNA counts; no whole-atlas dense conversion')
  if(is.null(rownames(counts))||is.null(colnames(counts))||anyDuplicated(rownames(counts))||anyDuplicated(colnames(counts))) stop('Counts axes require unique names')
  if(any(!is.finite(counts@x))||any(counts@x<0)||any(counts@x!=floor(counts@x)))stop('Not raw nonnegative integer RNA counts')
  total <- Matrix::colSums(counts)
  if(any(total<=0))stop('Zero RNA library; quarantine before annotation')
  x <- counts %*% Matrix::Diagonal(x=target_sum/total)
  x@x <- log1p(x@x)
  dimnames(x) <- dimnames(counts)
  x
}

display_umap_checked <- function(object,reduction='harmony',dims=1:20,
                                 name='UMAP_vNext',seed=20261005L) {
  emb <- SeuratObject::Embeddings(object,reduction=reduction)
  if(max(dims)>ncol(emb)||!all(is.finite(emb)))stop('Invalid reduction dimensions')
  if(!identical(rownames(emb),SeuratObject::Cells(object)))stop('Reduction cell axis differs from object')
  Seurat::RunUMAP(object,reduction=reduction,dims=dims,reduction.name=name,
                 seed.use=seed,n.neighbors=30L,min.dist=0.3,umap.method='uwot',
                 metric='cosine',return.model=TRUE)
}

write_reproducibility_record <- function(path,inputs,settings) {
  if(!requireNamespace('jsonlite',quietly=TRUE))stop('jsonlite required')
  jsonlite::write_json(list(input_files=normalizePath(inputs,mustWork=TRUE),
                           settings=settings,sessionInfo=capture.output(sessionInfo()),
                           note='UMAP is a display. Assess integration using high-dimensional embeddings and independent biological evidence.'),
                       path,auto_unbox=TRUE,pretty=TRUE)
}

export_RNA_layer_axes_checked <- function(object,path,source_file) {
  # Call on the exact object being exported, before matrix serialization.
  # Never derive dense/scale feature labels by VariableFeatures membership.
  counts<-get_RNA_counts_checked(object)
  if(any(!is.finite(counts@x))||any(counts@x<0)||any(counts@x!=floor(counts@x)))stop('RNA/counts does not contain raw nonnegative integers')
  scale<-if(utils::packageVersion('SeuratObject')>='5.0.0')
    SeuratObject::GetAssayData(object,assay='RNA',layer='scale.data') else
    SeuratObject::GetAssayData(object,assay='RNA',slot='scale.data')
  raw_axis<-rownames(counts);scaled_axis<-rownames(scale)
  stopifnot(!anyDuplicated(raw_axis),!anyDuplicated(scaled_axis))
  if(nrow(scale)&&!identical(colnames(scale),colnames(counts)))stop('Scaled RNA cell axis differs from raw RNA')
  if(nrow(scale)&&!all(scaled_axis%in%raw_axis))stop('Scaled feature names differ from raw RNA; explicitly reconcile source aliases before export')
  rna<-object[['RNA']]
  # The public getter maps Assay5's internal metadata rows to feature IDs;
  # its physical meta.data slot may use positional row names.
  md<-rna[[]]
  if(nrow(md)&&!is.null(rownames(md))) {
    raw_md<-md[match(raw_axis,rownames(md)),,drop=FALSE]
    scaled_md<-md[match(scaled_axis,rownames(md)),,drop=FALSE]
    if(!identical(rownames(raw_md),raw_axis)||!identical(rownames(scaled_md),scaled_axis))stop('Feature metadata does not align with matrix feature axes')
  }
  if(!requireNamespace('digest',quietly=TRUE))stop('digest required for source SHA256')
  jsonlite::write_json(list(schema_version='SCAID-RNA-source-axis-1.0',source_file=normalizePath(source_file,mustWork=TRUE),
                           source_file_sha256=digest::digest(file=source_file,algo='sha256'),
                           raw_feature_names=raw_axis,scaled_HVG_feature_names=scaled_axis,
                           cell_names=colnames(counts),source_feature_order_verified=TRUE,
                           export_association='This sidecar must be created from the exact object passed to the matrix exporter.',
                           raw_semantics='unintegrated nonnegative integer RNA counts',
                           scaled_HVG_semantics='RNA scale.data; not full log-normalized RNA'),
                       path,auto_unbox=TRUE,pretty=TRUE)
  invisible(list(raw_features=length(raw_axis),scaled_features=length(scaled_axis),cells=ncol(counts)))
}
