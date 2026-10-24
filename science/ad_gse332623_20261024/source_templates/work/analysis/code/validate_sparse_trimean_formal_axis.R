## Independent observed + two permutation means only; never infer networks.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
suppressPackageStartupMessages({library(CellChat);library(Matrix);library(hdf5r);library(future);library(jsonlite)})
work<-args[1];out<-args[2];dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
future::plan('sequential');if(requireNamespace('data.table',quietly=TRUE))data.table::setDTthreads(1L)
source(file.path(work,'pathway_and_cellchat/code/cellchat_exact_mean.R'))
source(file.path(work,'analysis/code/sparse_trimean_prototype.R'))
load_sparse_trimean_prototype(file.path(work,'analysis/code/sparse_trimean_prototype.cpp'),file.path(dirname(out),'pilot100_rcpp_cache'))
wait_memory<-function() {
  repeat {
    line<-grep('^MemAvailable:',readLines('/proc/meminfo'),value=TRUE)
    available<-as.numeric(strsplit(trimws(line),'[[:space:]]+')[[1]][2])*1024
    if(available>=320*1024^3)return(invisible(available))
    cat('WAIT_SYSTEM_MEMORY_RESERVE_320GiB\n');flush.console();Sys.sleep(10)
  }
}
wait_memory()
receipt_path<-file.path(work,'pathway_and_cellchat/performance_cache/normalized_DB_expression_cache_receipt.json')
r<-fromJSON(receipt_path);stopifnot(r$status=='validated_internal_cache_not_a_model_run',identical(as.integer(r$shape),c(581725L,1434L)))
f<-H5File$new(r$path,'r')
cells<-as.character(f[['obs/cell_id']][]);genes<-as.character(f[['gene/canonical_symbols']][])
labels<-as.character(f[['obs/MajorCellType']][]);keep<-as.logical(f[['obs/include_cellchat']][])
stopifnot(length(cells)==581725L,length(genes)==1434L,sum(keep)==528739L,!anyDuplicated(cells),!anyDuplicated(genes))
x<-f[['expression_csc/data']][];i<-as.integer(f[['expression_csc/indices']][]);p<-as.integer(f[['expression_csc/indptr']][])
f$close_all()
mc<-new('dgCMatrix',i=i,x=as.numeric(x),p=p,Dim=c(length(cells),length(genes)),Dimnames=list(cells,genes))
rm(x,i,p);gc(FALSE)
mc<-mc[keep,,drop=FALSE]
wait_memory();mat<-as.matrix(mc);rm(mc);gc(FALSE)
normalizer<-max(mat);stopifnot(is.finite(normalizer),normalizer>0)
mat<-mat/normalizer;gc(FALSE)
g<-factor(labels[keep],levels=sort(unique(labels[keep])))
stopifnot(nlevels(g)==8L,all(is.finite(mat)),all(mat>=0))
cat(toJSON(list(stage='formal_axis_ready',cells=nrow(mat),genes=ncol(mat),groups=levels(g),normalizer=normalizer),auto_unbox=TRUE),'\n');flush.console()
wait_memory();a<-proc.time()[['elapsed']];prepared<-scaid_prepare_positive_order(mat);prepare_seconds<-proc.time()[['elapsed']]-a
set.seed(2026L);permutation<-replicate(2L,sample.int(nrow(mat),size=nrow(mat)))
results<-list();saved<-list()
for(j in 0:2) {
  wait_memory();group<-if(j==0L)g else g[permutation[,j]]
  a<-proc.time()[['elapsed']];original<-cc_exact_group_aggregate(mat,group,CellChat::triMean);official_seconds<-proc.time()[['elapsed']]-a
  a<-proc.time()[['elapsed']];alternative<-cc_positive_group_aggregate(prepared,group,CellChat::triMean);sparse_seconds<-proc.time()[['elapsed']]-a
  stopifnot(identical(original,alternative))
  results[[j+1L]]<-list(grouping=if(j==0L)'Observed' else paste0('Seed2026_permutation_',j),
      gene_group_values=length(genes)*nlevels(g),bitwise_all_values_and_attributes_identical=TRUE,
      official_triMean_helper_seconds=official_seconds,sparse_positive_seconds=sparse_seconds)
  saved[[j+1L]]<-list(official=original,sparse=alternative)
  cat(toJSON(results[[j+1L]],auto_unbox=TRUE),'\n');flush.console()
}
saveRDS(list(results=saved,first_two_permutations=permutation,group_counts=table(g),normalizer=normalizer,genes=genes),paste0(out,'.rds'))
record<-list(status='validated_formal_axis_means_only_not_network',source_RNA_sha256=r$source_RNA_sha256,
 source_cache_sha256=r$sha256,cells=nrow(mat),genes=ncol(mat),groups=nlevels(g),group_counts=as.list(table(g)),
 prepare_seconds=prepare_seconds,comparisons=results,seed=2026L,permutations_checked=2L,
 scope='All 1434 genes and8groups on actual528739 CellChat records; normalized RNA unchanged, same max scaling and exact seed2026 first2 sample.int permutations; no communication model or running formal job modified',
 peak_memory_status=grep('^VmHWM:',readLines('/proc/self/status'),value=TRUE))
write_json(record,paste0(out,'.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('FORMAL_AXIS_THREE_GROUPINGS_BITWISE_VALIDATED\n')
