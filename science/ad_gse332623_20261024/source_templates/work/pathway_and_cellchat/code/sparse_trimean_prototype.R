## Independent prototype only; never sourced by either running formal job.
load_sparse_trimean_prototype<-function(cpp,cache_dir) {
  dir.create(cache_dir,recursive=TRUE,showWarnings=FALSE)
  Rcpp::sourceCpp(cpp,cacheDir=cache_dir,showOutput=FALSE,verbose=FALSE)
}
cc_positive_group_aggregate<-function(prepared,group,FUN) {
  stopifnot(identical(FUN,CellChat::triMean),is.factor(group),
            !anyNA(group),nlevels(group)==length(unique(group)))
  labels<-levels(group)
  q<-scaid_positive_group_quartiles(prepared,as.integer(group),length(labels))
  means<-apply(q,c(2L,3L),mean)
  answer<-data.frame(Group.1=factor(labels,levels=labels),t(means),
                     check.names=FALSE,row.names=seq_along(labels))
  colnames(answer)<-c('Group.1',prepared$gene_names)
  answer
}
make_positive_order_CellChat<-function() {
  fun<-make_exact_mean_serial_LR_CellChat()
  changed_prep<-0L;changed_aggregate<-0L
  replace<-function(expr) {
    if(!is.call(expr))return(expr)
    if(identical(expr,quote(data.use.t<-t(data.use)))) {
      changed_prep<<-changed_prep+1L
      return(quote(data.use.positive<-.cc_prepare_positive_order(t(data.use))))
    }
    if(identical(expr[[1]],as.name('.cc_exact_group_aggregate'))) {
      changed_aggregate<<-changed_aggregate+1L
      expr[[1]]<-as.name('.cc_positive_group_aggregate')
      stopifnot(identical(expr[[2]],as.name('data.use.t')))
      expr[[2]]<-as.name('data.use.positive')
    }
    for(i in seq_along(expr))if(i>1L&&is.call(expr[[i]]))expr[[i]]<-replace(expr[[i]])
    expr
  }
  body(fun)<-replace(body(fun))
  stopifnot(changed_prep==1L,changed_aggregate==2L)
  assign('.cc_prepare_positive_order',scaid_prepare_positive_order,envir=environment(fun))
  assign('.cc_positive_group_aggregate',cc_positive_group_aggregate,envir=environment(fun))
  fun
}
