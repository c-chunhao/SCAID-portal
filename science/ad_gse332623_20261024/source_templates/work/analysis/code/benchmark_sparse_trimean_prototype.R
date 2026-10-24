## Independent prototype verification only. Does not start a formal candidate.
args<-commandArgs(TRUE);stopifnot(length(args)==4L)
suppressPackageStartupMessages({library(CellChat);library(future);library(jsonlite)})
future::plan('sequential')
work<-args[1];out<-args[4];dir.create(dirname(out),recursive=TRUE,showWarnings=FALSE)
source(file.path(work,'pathway_and_cellchat/code/cellchat_exact_mean.R'))
source(file.path(work,'analysis/code/sparse_trimean_prototype.R'))
load_sparse_trimean_prototype(file.path(work,'analysis/code/sparse_trimean_prototype.cpp'),paste0(out,'_rcpp_cache'))
set.seed(8037)
edge_cases<-0L;edge_outputs<-0L
for(n in c(1:80,101L,256L,347L,1001L)) {
  ng<-min(8L,n);g<-factor(rep(paste0('Group',seq_len(ng)),length.out=n),levels=paste0('Group',seq_len(ng)))
  m<-matrix(runif(n*16),nrow=n);m[,1]<-0;m[,2]<-.1;m[,3]<-1
  for(j in 4:16)m[runif(n)<((j-4)/13),j]<-0
  m[,4]<-sample(c(0,0,0,1e-200,1e-16,.1,.3,1),n,replace=TRUE)
  colnames(m)<-paste0('Gene',1:16)
  p<-scaid_prepare_positive_order(m)
  standard<-aggregate(m,list(g),FUN=CellChat::triMean)
  alternative<-cc_positive_group_aggregate(p,g,CellChat::triMean)
  stopifnot(identical(standard,alternative))
  edge_cases<-edge_cases+1L;edge_outputs<-edge_outputs+16L*ng
}
set.seed(2026)
n<-347L;g<-factor(rep(LETTERS[1:8],length.out=n),levels=LETTERS[1:8])
m<-matrix(sample(c(0,0,0,0,.1,.1,.2,.3,1,9,1e-200),n*17,replace=TRUE),nrow=n)
m[,1]<-0;m[,2]<-1;colnames(m)<-c('ALL_ZERO','ALL_TIED',paste0('Gene',3:17));p<-scaid_prepare_positive_order(m)
perm<-replicate(100L,sample.int(n,n))
for(k in seq_len(ncol(perm))) {
  b<-g[perm[,k]]
  stopifnot(identical(aggregate(m,list(b),FUN=CellChat::triMean),cc_positive_group_aggregate(p,b,CellChat::triMean)))
}
bad<-m;bad[1,1]<-NA_real_
stopifnot(inherits(try(scaid_prepare_positive_order(bad),silent=TRUE),'try-error'))
bad[1,1]<--1
stopifnot(inherits(try(scaid_prepare_positive_order(bad),silent=TRUE),'try-error'))
cat('SPARSE_EDGE_TESTS_PASSED\n');flush.console()
object<-readRDS(args[2]);reference<-readRDS(args[3]);stopifnot(ncol(object@data)==17536L,nlevels(object@idents)==7L)
fast<-make_positive_order_CellChat()
phase<-new.env(parent=emptyenv());phase$prepare_seconds<-0;phase$aggregate_seconds<-0;phase$aggregate_calls<-0L
env<-environment(fast);oldprep<-get('.cc_prepare_positive_order',envir=env);oldaggregate<-get('.cc_positive_group_aggregate',envir=env)
assign('.cc_prepare_positive_order',function(x){a<-proc.time()[['elapsed']];v<-oldprep(x);phase$prepare_seconds<-phase$prepare_seconds+proc.time()[['elapsed']]-a;v},envir=env)
assign('.cc_positive_group_aggregate',function(p,g,FUN){a<-proc.time()[['elapsed']];v<-oldaggregate(p,g,FUN);phase$aggregate_seconds<-phase$aggregate_seconds+proc.time()[['elapsed']]-a;phase$aggregate_calls<-phase$aggregate_calls+1L;v},envir=env)
elapsed<-system.time({result<-fast(object,type='triMean',raw.use=TRUE,population.size=FALSE,distance.use=FALSE,nboot=100L,seed.use=2026L)})[['elapsed']]
stopifnot(identical(reference$official_probability,result@net$prob),identical(reference$official_pvalue,result@net$pval),phase$aggregate_calls==101L)
saveRDS(list(prob=result@net$prob,pval=result@net$pval),paste0(out,'.rds'))
receipt<-list(status='prototype_bitwise_equivalent_not_formal_release',scope='Independent nonnegative finite sparse-positive quartile prototype; both running formal jobs and all source H5 files untouched',
 cells=ncol(object@data),groups=nlevels(object@idents),bootstrap=100L,seed=2026L,
 sparse_positive_seconds=elapsed,prepare_seconds=phase$prepare_seconds,aggregate_seconds=phase$aggregate_seconds,aggregate_calls=phase$aggregate_calls,
 probability_max_absolute_error=max(abs(reference$official_probability-result@net$prob)),pvalue_max_absolute_error=max(abs(reference$official_pvalue-result@net$pval)),
 all_probability_pvalue_and_dimnames_bitwise_identical=TRUE,
 edge_test=list(cases=edge_cases,group_gene_outputs=edge_outputs,permutations=100L,additional_group_gene_outputs=13600L,n1_to_80_and_large_odd_even=TRUE,all_zero_and_ties=TRUE,very_small_and_arbitrary_double=TRUE,NA_and_negative_explicitly_rejected=TRUE),
 arithmetic='collapse2.1.6 type7 a+gamma*(b-a), with explicit double operation order; official R mean combines Q25,Q50,Q50,Q75 unchanged',
 packages=as.list(sapply(c('CellChat','collapse','Rcpp'),function(p)as.character(packageVersion(p)))))
write_json(receipt,paste0(out,'.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('SPARSE_FULL_100BOOT_PILOT_PASSED\n')
