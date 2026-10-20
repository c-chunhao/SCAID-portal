#!/usr/bin/env Rscript
# Export the existing processed Seurat object without filtering or relabeling it.
suppressPackageStartupMessages({library(Seurat); library(Matrix); library(rhdf5); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
if(!(length(args)%in%c(1L,2L,3L))) stop('Usage: export_new_aid_h5ad.R object_index [work_root] [recover_raw_counts]')
object_index <- as.integer(args[[1]])
work <- if(length(args)>=2L) args[[2]] else Sys.getenv('SCAID_SCIENCE_WORK_ROOT','work')
dir.create(work,recursive=TRUE,showWarnings=FALSE)
work <- normalizePath(work,mustWork=TRUE)
plan_path <- Sys.getenv('SCAID_RELEASE_PLAN','examples/release_plan.example.json')
plan <- fromJSON(plan_path,simplifyVector=FALSE)
if(object_index<0L||object_index>=length(plan$objects)) stop('Object index outside release plan')
record <- plan$objects[[object_index+1L]]
object_id <- record$object_id
recover_raw <- length(args)==3L&&args[[3]]=='recover_raw_counts'
if(recover_raw&&object_id!='T1D/GSE279086/Kidney') stop('Special recovery policy is restricted to independently audited GSE279086')
stem <- gsub('/','_',object_id,fixed=TRUE)
input_stem <- gsub('/','-',object_id,fixed=TRUE)
input_root <- Sys.getenv('SCAID_SEURAT_INPUT_ROOT',file.path(Sys.getenv('SCAID_SCIENCE_INPUT_ROOT','inputs'),'00.data'))
input <- file.path(input_root,object_id,paste0(input_stem,'.rds'))
outroot <- Sys.getenv('SCAID_EXPORT_ROOT',file.path(work,'exports'))
reportdir <- file.path(work,'matrices','objects',stem)
dir.create(reportdir,recursive=TRUE,showWarnings=FALSE)
report_path <- file.path(reportdir,'export_report.json')
started <- format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z')
message('START ',object_id,' ',started)

attr_set <- function(file,path,attrs) {
  fid <- H5Fopen(file); obj <- H5Oopen(fid,path)
  on.exit({H5Oclose(obj);H5Fclose(fid)})
  for(n in names(attrs)) h5writeAttribute(if(is.logical(attrs[[n]])) as.integer(attrs[[n]]) else attrs[[n]],obj,n,asScalar=length(attrs[[n]])==1L&&!(n%in%c('column-order','shape')),
                                       encoding='UTF-8',variableLengthString=is.character(attrs[[n]]))
}
encoded_group <- function(file,path,type='dict',version='0.1.0') {
  if(path!='/') stopifnot(h5createGroup(file,path))
  attr_set(file,path,list('encoding-type'=type,'encoding-version'=version))
}
array_write <- function(file,path,x) {
  d <- dim(x)
  chunk <- if(is.null(d)) max(1L,min(length(x),1000000L)) else if(length(d)==2L) c(d[[1]],min(d[[2]],2000L)) else d
  if(is.character(x)) {
    h5write(enc2utf8(x),file,path,encoding='UTF-8',level=1)
    attr_set(file,path,list('encoding-type'='string-array','encoding-version'='0.2.0'))
  } else {
    h5createDataset(file,path,dims=if(is.null(d)) length(x) else d,
                    storage.mode=if(is.integer(x)) 'integer' else 'double',chunk=chunk,level=1)
    h5write(x,file,path)
    attr_set(file,path,list('encoding-type'='array','encoding-version'='0.2.0'))
  }
}
scalar_write <- function(file,path,x) {
  h5write(x,file,path)
  # Scalar strings are written as a one-element array by rhdf5. Keep the standard array encoding.
  attr_set(file,path,list('encoding-type'=if(is.character(x)) 'string-array' else 'array','encoding-version'='0.2.0'))
}
frame_write <- function(file,path,frame,index) {
  stopifnot(nrow(frame)==length(index),!anyDuplicated(index),!anyNA(index),!any(grepl('/',names(frame),fixed=TRUE)),!('_index'%in%names(frame)))
  encoded_group(file,path,'dataframe','0.2.0')
  attr_set(file,path,list('_index'='_index','column-order'=as.character(names(frame))))
  array_write(file,paste0(path,'/_index'),index)
  for(n in names(frame)) {
    x <- frame[[n]]; p <- paste0(path,'/',n)
    if(is.character(x)||is.factor(x)||is.logical(x)) {
      lev <- if(is.factor(x)) levels(x) else unique(as.character(x[!is.na(x)]))
      # A nonempty category table also supports entirely missing source columns.
      if(length(lev)==0L) lev <- '__no_observed_value__'
      codes <- match(as.character(x),lev)-1L; codes[is.na(x)] <- -1L
      encoded_group(file,p,'categorical','0.2.0')
      attr_set(file,p,list(ordered=FALSE))
      array_write(file,paste0(p,'/codes'),as.integer(codes))
      array_write(file,paste0(p,'/categories'),as.character(lev))
    } else if(is.integer(x)||is.numeric(x)) {
      # Numeric NA is written as NaN, avoiding a sentinel mistaken for a count.
      x <- as.numeric(x); x[is.na(x)] <- NaN
      array_write(file,p,x)
    } else stop('Unsupported source metadata class ',n,': ',paste(class(x),collapse=','))
  }
}
array_equal <- function(a,b) isTRUE(all.equal(as.vector(a),as.vector(b),check.attributes=FALSE,tolerance=0))
vector_verify <- function(file,path,x,chunk=2000000L) {
  n <- length(x)
  if(n==0L) return(TRUE)
  for(s in seq.int(1L,n,by=chunk)) {
    k <- min(chunk,n-s+1L)
    y <- h5read(file,path,start=s,count=k)
    if(!array_equal(x[seq.int(s,length.out=k)],y)) stop('Full readback differs: ',path,' at ',s)
  }
  TRUE
}
sparse_write_verify <- function(file,path,m) {
  if(!inherits(m,'dgCMatrix')) m <- as(m,'dgCMatrix')
  if(any(!is.finite(m@x))) stop('Non-finite sparse values: ',path)
  encoded_group(file,path,'csr_matrix','0.1.0')
  attr_set(file,path,list(shape=as.integer(rev(dim(m)))))
  array_write(file,paste0(path,'/data'),m@x)
  array_write(file,paste0(path,'/indices'),m@i)
  array_write(file,paste0(path,'/indptr'),m@p)
  # Seurat is genes x cells CSC. These same components encode cells x genes CSR.
  for(n in c('data','indices','indptr')) vector_verify(file,paste0(path,'/',n),switch(n,data=m@x,indices=m@i,indptr=m@p))
  list(source_class=class(m)[1],shape=as.integer(rev(dim(m))),nnz=length(m@x),
       minimum=if(length(m@x)) min(0,m@x) else 0,maximum=if(length(m@x)) max(0,m@x) else 0,
       finite=TRUE,full_component_readback=TRUE)
}
dense_write_verify <- function(file,path,m) {
  stopifnot(is.matrix(m))
  if(any(!is.finite(m))) stop('Non-finite dense matrix: ',path)
  # rhdf5 reverses the two dimensions when exposed through HDF5; original
  # features x cells therefore becomes the required AnnData cells x features.
  array_write(file,path,m)
  for(s in seq.int(1L,ncol(m),by=2000L)) {
    k <- min(2000L,ncol(m)-s+1L)
    y <- h5read(file,path,start=c(1L,s),count=c(nrow(m),k))
    if(!array_equal(m[,seq.int(s,length.out=k),drop=FALSE],y)) stop('Dense readback differs: ',path,' at ',s)
  }
  list(source_class=class(m)[1],shape=as.integer(rev(dim(m))),minimum=min(m),maximum=max(m),finite=TRUE,full_component_readback=TRUE)
}
matrix_write_verify <- function(file,path,m) {
  z <- if(inherits(m,'sparseMatrix')) sparse_write_verify(file,path,m) else dense_write_verify(file,path,as.matrix(m))
  rr <- unique(as.integer(round(seq.int(1L,nrow(m),length.out=min(17L,nrow(m))))))
  cc <- unique(as.integer(round(seq.int(1L,ncol(m),length.out=min(19L,ncol(m))))))
  ij <- expand.grid(feature=rr,cell=cc)
  z$source_value_probes <- list(cell_index_zero=as.integer(ij$cell-1L),feature_index_zero=as.integer(ij$feature-1L),values=as.numeric(m[cbind(ij$feature,ij$cell)]))
  z
}
layer_get <- function(o,a,layer) {
  z <- LayerData(o,assay=a,layer=layer,fast=FALSE)
  if(!identical(colnames(z),colnames(o))) stop(a,'/',layer,' cell axis differs from source object')
  if(anyDuplicated(rownames(z))||anyNA(rownames(z))) stop(a,'/',layer,' duplicate or missing features')
  z
}
summary_column <- function(x) {
  z <- as.character(x); valid <- !is.na(z)&nzchar(trimws(z))&!(tolower(trimws(z))%in%c('na','nan','null','none','unknown'))
  list(observed=sum(valid),missing=sum(!valid),unique_observed=length(unique(z[valid])),
       levels=as.list(sort(table(z[valid]),decreasing=TRUE)))
}
result <- tryCatch({
  o <- readRDS(input)
  stopifnot(inherits(o,'Seurat'),ncol(o)==as.integer(record$catalog$cell_count),
            !anyDuplicated(colnames(o)),identical(rownames(o@meta.data),colnames(o)))
  if(!all(c('RNA','AUCell','UCell','singscore')%in%Assays(o))) stop('One or more required assays absent')
  meta <- o@meta.data
  # Original source columns remain unchanged. Portal identifiers are additive.
  meta$scaid_object_id <- object_id
  meta$scaid_condition <- record$catalog$disease_name
  meta$scaid_source <- record$catalog$dataset_id
  meta$scaid_tissue <- record$catalog$tissue
  stopifnot(!anyDuplicated(names(meta)))
  emb <- list()
  for(n in c('pca','harmony','umap')) if(n%in%Reductions(o)) {
    e <- Embeddings(o,n)
    if(!identical(rownames(e),colnames(o))||any(!is.finite(e))) stop('Invalid reduction ',n)
    emb[[paste0('X_',n)]] <- e
  }
  md <- lapply(meta[,intersect(c('Sample','GSM','Group','Disease_detail','Disease_subtype','MajorCellType','Immune_All_High','doublet_status','Tissue','Dataset'),names(meta)),drop=FALSE],summary_column)
  samples <- if('Sample'%in%names(meta)) as.data.frame(table(meta$Sample,useNA='ifany'),stringsAsFactors=FALSE) else data.frame()
  if(nrow(samples)) {names(samples)<-c('original_sample_label','cell_count');write.table(samples,file.path(reportdir,'sample_labels.tsv'),sep='\t',row.names=FALSE,quote=TRUE,na='NA')}
  cells <- data.frame(cell_id=colnames(o),meta,check.names=FALSE)
  gz <- gzfile(file.path(reportdir,'cell_metadata.tsv.gz'),'wt');write.table(cells,gz,sep='\t',row.names=FALSE,quote=TRUE,na='__SCAID_SOURCE_MISSING__');close(gz);rm(cells)
  matrices <- list()
  for(a in c('RNA','AUCell','singscore','UCell')) {
    dir.create(file.path(outroot,a),recursive=TRUE,showWarnings=FALSE)
    final <- file.path(outroot,a,paste0(stem,'_',a,'.h5ad'))
    temp <- paste0(final,'.tmp')
    if(file.exists(final)) stop('Refusing to replace existing validated output: ',final)
    if(file.exists(temp)) unlink(temp)
    message('EXPORT ',object_id,' ',a)
    stopifnot(h5createFile(temp))
    encoded_group(temp,'/','anndata','0.1.0')
    m <- layer_get(o,a,'data')
    features <- rownames(m)
    frame_write(temp,'obs',meta,colnames(o))
    feature_frame <- data.frame(feature_name=features,stringsAsFactors=FALSE)
    if(a!='RNA') {feature_frame$kegg_id<-sub('^(hsa[0-9]+).*','\\1',features);feature_frame$score_method<-a}
    frame_write(temp,'var',feature_frame,features)
    for(n in c('obsm','varm','obsp','varp','layers','uns')) encoded_group(temp,n)
    for(n in names(emb)) array_write(temp,paste0('obsm/',n),t(emb[[n]]))
    encoded_group(temp,'uns/scaid')
    provenance <- list(object_id=object_id,source_accession=record$catalog$dataset_id,
      source_rds=input,source_assay=a,source_layer='data',
      X_representation=if(a=='RNA') 'RNA data layer preserved from the processed Seurat object' else paste(a,'KEGG pathway scores; original unscaled data layer'),
      counts_representation=if(a=='RNA'&&recover_raw) 'Integer counts recovered by exact cell/feature matching from source raw matrices; original fractional processed matrix is preserved separately as source_processed_counts' else if(a=='RNA') 'Original RNA counts layer, preserved without rounding or filtering' else 'Not applicable: a pathway score matrix',
      sample_labels='Source sample labels are preserved; sample labels are not authenticated donor IDs',
      processing='Format conversion only; no additional filtering, annotation, normalization or batch correction',
      exported_at=started)
    for(n in names(provenance)) scalar_write(temp,paste0('uns/scaid/',n),provenance[[n]])
    v <- list(X=matrix_write_verify(temp,'X',m))
    if(a=='RNA') {
      counts <- layer_get(o,'RNA','counts')
      if(!identical(rownames(counts),features)) stop('RNA counts feature axis differs from data')
      if(any(counts@x<0)||any(!is.finite(counts@x))) stop('RNA counts are not finite and nonnegative')
      if(any(counts@x!=round(counts@x))) {
        if(!recover_raw) stop('RNA counts are not finite nonnegative integers')
        v$source_processed_counts <- matrix_write_verify(temp,'layers/source_processed_counts',counts)
        scalar_write(temp,'uns/scaid/processed_counts_note','The source RNA counts layer contains fractional processed expression values; preserved as layers/source_processed_counts, not represented as raw counts. X remains the original normalized source data layer.')
        scalar_write(temp,'uns/scaid/nCount_RNA_semantics','Original nCount_RNA and nFeature_RNA describe source_processed_counts; separately recovered integer raw counts have additional *_source_raw metadata.')
      } else {
        v$counts <- matrix_write_verify(temp,'layers/counts',counts)
        v$counts$nonnegative_integer_values <- TRUE
      }
      rm(counts)
    }
    stopifnot(identical(as.vector(h5read(temp,'obs/_index')),colnames(o)),identical(as.vector(h5read(temp,'var/_index')),features))
    for(n in names(emb)) stopifnot(array_equal(h5read(temp,paste0('obsm/',n)),t(emb[[n]])))
    H5close()
    matrices[[a]] <- list(method=a,path=final,temporary_path=temp,shape=as.integer(c(ncol(o),length(features))),
      feature_count=length(features),cell_count=ncol(o),source_layers=Layers(o[[a]]),
      full_source_component_readback=TRUE,validated=FALSE,validation_status='awaiting_independent_anndata_readback',matrix_checks=v)
    rm(m);gc(verbose=FALSE)
  }
  list(status=if(recover_raw) 'exported_awaiting_raw_count_recovery' else 'exported_awaiting_python_validation',object_index=object_index,object_id=object_id,
       source_rds=input,source_bytes=unname(file.info(input)$size),source_mtime=as.character(file.info(input)$mtime),
       started=started,finished=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'),catalog_cell_count=record$catalog$cell_count,
       actual_cell_count=ncol(o),actual_rna_genes=nrow(o[['RNA']]),metadata_columns=names(meta),metadata_summary=md,
       reductions=lapply(emb,dim),matrices=matrices,
       versions=list(R=as.character(getRversion()),Seurat=as.character(packageVersion('Seurat')),SeuratObject=as.character(packageVersion('SeuratObject')),rhdf5=as.character(packageVersion('rhdf5'))))
},error=function(e) {message('ERROR ',object_id,' ',conditionMessage(e));list(status='failed',object_index=object_index,object_id=object_id,source_rds=input,error=conditionMessage(e),started=started,finished=format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'))})
write_json(result,paste0(report_path,'.tmp'),pretty=TRUE,auto_unbox=TRUE,na='null',null='null',digits=NA)
stopifnot(file.rename(paste0(report_path,'.tmp'),report_path))
H5close()
if(result$status=='failed') quit(status=1L)
message('EXPORTED ',object_id)
