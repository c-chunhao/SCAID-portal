local({
## One-library, checkpointed counts-only ingest. No QC, normalization or inferred clinical groups.
## Usage: Rscript code/GSE332623_build_corrected_input.R <library_id> [output_root]
args <- commandArgs(TRUE)
if(length(args)<1L || length(args)>2L) stop('Usage: <library_id> [output_root]')
library_id<-args[1L]
output_root<-if(length(args)==2L)args[2L]else'@PROJECT@/00.corrected_data/GSE332623'
if(!grepl('^GSE332623_TotalseqB_NEA[1-7]_ADT[1-4]$',library_id))stop('Invalid library ID')
if(basename(normalizePath(output_root,mustWork=FALSE))!='GSE332623')stop('Output must be a dedicated GSE332623 directory')
source('@PROJECT@/code/v2/GSE332623.multimodal.R')
suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(jsonlite)})
input_dir<-'@RAW_INPUT@'
inv_path<-'@PROJECT@/code/v2/00.Reference/SampleInventory.tsv'
inv<-read.delim(inv_path,sep='\t',quote='',colClasses='character',check.names=FALSE)
dir.create(file.path(output_root,'libraries'),recursive=TRUE,showWarnings=FALSE)
dir.create(file.path(output_root,'receipts'),recursive=TRUE,showWarnings=FALSE)
output<-file.path(output_root,'libraries',paste0(library_id,'.rds'))
receipt_path<-file.path(output_root,'receipts',paste0(library_id,'.json'))
sha256 <- function(path) {
 result<-system2('sha256sum',shQuote(path),stdout=TRUE,stderr=TRUE)
 if(!is.null(attr(result,'status'))||length(result)!=1L)stop('Hash failed')
 substr(result,1L,64L)
}
helper_path<-'@PROJECT@/code/v2/GSE332623.multimodal.R'
helper_sha<-sha256(helper_path)
if(file.exists(receipt_path) && file.exists(output)) {
 previous<-read_json(receipt_path,simplifyVector=TRUE)
 if(identical(previous$status,'validated') && identical(previous$reader_sha256,helper_sha) &&
    identical(as.numeric(file.info(output)$size),as.numeric(previous$output_size_bytes)) &&
    identical(sha256(output),previous$output_sha256)){
  paths<-unlist(previous$source_files); hashes<-unlist(previous$source_sha256)
  if(all(vapply(seq_along(paths),function(i)identical(sha256(paths[i]),unname(hashes[i])),logical(1)))){
   cat('SKIP_VALIDATED',library_id,'\n');quit(save='no',status=0L)
  }
 }
 stop('Existing input/receipt differs: inspect it before replacing ',library_id)
}
if(file.exists(output)||file.exists(receipt_path))stop('Incomplete previous checkpoint exists; inspect before replacing')
started<-Sys.time();cat('READ_START',library_id,format(started),'\n')
obj<-ReadGSE332623Library(input_dir,library_id,inv)
source_input<-obj@misc$source_input
source_files<-unlist(source_input$paths)
source_sha<-vapply(source_files,sha256,character(1))
expected_cells<-colnames(obj)
expected_rna<-rownames(obj[['RNA']]);expected_adt<-rownames(obj[['ADT']])
expected_rna_x<-LayerData(obj,assay='RNA',layer='counts');expected_adt_x<-LayerData(obj,assay='ADT',layer='counts')
# A numerical fingerprint supplements exact feature/barcode axes and NNZ preservation.
fingerprint<-function(m)list(nnzero=length(m@x),sum=sum(m@x),sum_squared=sum(m@x*m@x),indptr_sum=sum(as.numeric(m@p)),index_sum=sum(as.numeric(m@i)))
expected_fingerprints<-list(RNA=fingerprint(expected_rna_x),ADT=fingerprint(expected_adt_x))
rm(expected_rna_x,expected_adt_x);invisible(gc())
tmp<-paste0(output,'.tmp.',Sys.getpid());con<-gzfile(tmp,'wb',compression=1L)
saveRDS(obj,con);close(con)
rm(obj);invisible(gc())
check<-readRDS(tmp)
stopifnot(identical(colnames(check),expected_cells),
          identical(rownames(check[['RNA']]),expected_rna),
          identical(rownames(check[['ADT']]),expected_adt),
          identical(colnames(check[['RNA']]),colnames(check[['ADT']])),
          all(check$Group=='Unresolved'),all(check$source_condition=='Unresolved'),
          all(is.na(check$donor_id)),length(check@reductions)==0L,
          setequal(Assays(check),c('RNA','ADT')),
          all(check[['RNA']][[]]$feature_type=='Gene Expression'),
          all(check[['ADT']][[]]$feature_type=='Antibody Capture'))
for(a in c('RNA','ADT')){
 stopifnot(identical(Layers(check[[a]]),'counts'))
 m<-LayerData(check,assay=a,layer='counts')
 stopifnot(identical(fingerprint(m),expected_fingerprints[[a]]),all(is.finite(m@x)),
           all(m@x>=0),all(m@x==round(m@x)))
}
if(!all(vapply(seq_along(source_files),function(i)identical(sha256(source_files[i]),unname(source_sha[i])),logical(1))))stop('Source files changed while reading')
output_sha<-sha256(tmp)
if(!file.rename(tmp,output))stop('Atomic checkpoint move failed')
receipt<-list(status='validated',stage='Unfiltered counts-only multimodal input; not analyzed AD cases',
 library_id=library_id,n_cells=ncol(check),rna_features=nrow(check[['RNA']]),adt_features=nrow(check[['ADT']]),
 source_files=as.list(source_files),source_sha256=as.list(source_sha),source_nnzero=source_input$source_nnzero,
 output_path=output,output_size_bytes=file.info(output)$size,output_sha256=output_sha,reader_sha256=helper_sha,
 source_gsm=unname(unique(check$source_gsm)),source_pool=unname(unique(check$source_pool)),
 source_condition='Unresolved',Group='Unresolved',donor_mapping='not available',
 retained_all_source_barcodes=as.numeric(ncol(check)) == as.numeric(source_input$source_barcode_count),
 fingerprints=expected_fingerprints,started_utc=format(started,tz='UTC',usetz=TRUE),finished_utc=format(Sys.time(),tz='UTC',usetz=TRUE),
 elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),
 versions=list(R=as.character(getRversion()),SeuratObject=as.character(packageVersion('SeuratObject')),Matrix=as.character(packageVersion('Matrix'))))
receipt_tmp<-paste0(receipt_path,'.tmp.',Sys.getpid());write_json(receipt,receipt_tmp,pretty=TRUE,auto_unbox=TRUE)
if(!file.rename(receipt_tmp,receipt_path))stop('Atomic receipt move failed')
cat('VALIDATED',library_id,ncol(check),nrow(check[['RNA']]),nrow(check[['ADT']]),'\n')

})
