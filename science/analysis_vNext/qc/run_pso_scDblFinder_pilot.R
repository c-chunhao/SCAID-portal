#!/usr/bin/env Rscript
args<-commandArgs(trailingOnly=TRUE)
if (!length(args) || identical(args[[1]],'--help')) {
  cat('Usage: Rscript run_pso_scDblFinder_pilot.R INPUT_ROOT [CALLING_STAGE] [ADAPTER_PATH]\n')
  cat('Frozen six-capture GSE162183 pilot replay; not a donor-assignment or generic-cohort pipeline.\n')
  quit(status=if (!length(args)) 2 else 0)
}
suppressPackageStartupMessages({
  library(Matrix);library(SingleCellExperiment);library(scDblFinder)
  library(BiocParallel);library(jsonlite);library(digest)
})
args<-commandArgs(trailingOnly=TRUE)
root<-args[[1]]
dir.create(file.path(root,'qc'),showWarnings=FALSE,recursive=TRUE)
dir.create(file.path(root,'provenance'),showWarnings=FALSE,recursive=TRUE)
calling_stage<-if(length(args)>1)args[[2]]else'before_distribution_tail_QC_primary'
file_arg<-grep('^--file=',commandArgs(),value=TRUE)[1]
script_dir<-dirname(normalizePath(sub('^--file=','',file_arg)))
adapter<-if(length(args)>2)args[[3]]else file.path(script_dir,'xgboost_modern_adapter.R')
stopifnot(as.character(packageVersion('scDblFinder'))=='1.16.0',as.character(packageVersion('xgboost'))=='3.2.1.1')
source(adapter)
channels<-c('Ctrl1_10x','Ctrl2_10x','Ctrl3_10x','Psor1_10x','Psor2_10x','Psor3_10x')
out<-list();summaries<-list();input_manifests<-list()
for(i in seq_along(channels)){
  ch<-channels[[i]];seed<-162183L+i
  input<-file.path(root,'doublets_input',paste0(ch,'.mtx'))
  ids_path<-file.path(root,'doublets_input',paste0(ch,'_barcodes.tsv'))
  ids<-read.delim(ids_path,stringsAsFactors=FALSE)$cell_id
  mat<-as(readMM(input),'CsparseMatrix')
  stopifnot(length(ids)==ncol(mat),!anyNA(ids),all(nzchar(ids)),!anyDuplicated(ids),all(is.finite(mat@x)),all(mat@x>=0),all(mat@x==round(mat@x)))
  colnames(mat)<-ids;rownames(mat)<-paste0('source_gene_',seq_len(nrow(mat)))
  sce<-SingleCellExperiment(assays=list(counts=mat))
  .adapter_audit$channel<-ch;.adapter_audit$probe<-NULL;.adapter_audit$previous<-NULL
  set.seed(seed)
  cat(ch,'CALLING',ncol(mat),'per-capture calling stage',calling_stage,'seed',seed,'\n');flush.console()
  sce<-scDblFinder(sce,clusters=FALSE,BPPARAM=SerialParam(RNGseed=seed))
  if(length(.adapter_audit$errors))stop(paste('Fatal classifier training error',paste(.adapter_audit$errors,collapse=' | ')))
  training<-do.call(rbind,.adapter_audit$calls);a<-training[training$channel==ch,]
  stopifnot(nrow(a)==3L,all(a$CV_executed),all(a$DMatrix_labels_exact),all(a$actual_trained_rounds>0L),all(a$probe_mean_change_from_previous[-1]>0))
  tab<-data.frame(cell_id=colnames(sce),capture_channel=ch,scDblFinder_score=sce$scDblFinder.score,
                  scDblFinder_class=as.character(sce$scDblFinder.class),seed=seed,stringsAsFactors=FALSE)
  stopifnot(nrow(tab)==length(ids),identical(tab$cell_id,ids),all(is.finite(tab$scDblFinder_score)),all(tab$scDblFinder_score>=0 & tab$scDblFinder_score<=1),all(tab$scDblFinder_class %in% c('singlet','doublet')))
  write.table(tab,file.path(root,'qc',paste0(ch,'_scDblFinder.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
  write.table(training,file.path(root,'provenance','xgboost_actual_training_audit.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
  out[[ch]]<-tab;summaries[[ch]]<-data.frame(capture_channel=ch,tested_cells=nrow(tab),doublets=sum(tab$scDblFinder_class=='doublet'),singlets=sum(tab$scDblFinder_class=='singlet'),seed=seed)
  input_manifests[[ch]]<-list(matrix_SHA256=digest(file=input,algo='sha256'),barcodes_SHA256=digest(file=ids_path,algo='sha256'),seed=seed,
                  scDblFinder=as.character(packageVersion('scDblFinder')),xgboost=as.character(packageVersion('xgboost')),
                  expected_rate='unmodified scDblFinder default from source processed-capture size, not fixed6%; loading rate remains unknown',original_counts_unchanged=TRUE)
  cat(ch,'DONE',summaries[[ch]]$doublets,'doublets',summaries[[ch]]$singlets,'singlets\n');flush.console()
  rm(mat,sce);gc()
}
write.table(do.call(rbind,out),file.path(root,'qc','scDblFinder_all_channels.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,summaries),file.path(root,'qc','scDblFinder_summary.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
sink(file.path(root,'provenance','R_sessionInfo.txt'));print(sessionInfo());sink()
write_json(list(inputs=input_manifests,adapter_path=adapter,adapter_SHA256=digest(file=adapter,algo='sha256'),
           strict_handler_count=.strict_handler_count,classifier_errors=.adapter_audit$errors,all_classifier_errors_fatal=TRUE,
           physical_capture_partition_verified=TRUE,calling_stage=calling_stage,source_preprocessing_history_not_replaced=TRUE,
           flags_only_no_release_mutation=TRUE),file.path(root,'provenance','doublet_execution_manifest.json'),auto_unbox=TRUE,pretty=TRUE)
cat('COMPLETE',sum(vapply(out,nrow,integer(1))),'source barcodes\n')
