## Standby prototype: identical installed CellChat::triMean and group order;
## avoid repeated whole-matrix transpose/data.frame construction per bootstrap.
cc_exact_group_aggregate <- function(transposed,group,FUN) {
  labels<-levels(group)
  value<-vapply(labels,function(label) {
    idx<-which(group==label)
    vapply(seq_len(ncol(transposed)),function(g)FUN(transposed[idx,g]),numeric(1))
  },numeric(ncol(transposed)))
  ans<-data.frame(Group.1=factor(labels,levels=labels),t(value),check.names=FALSE,row.names=seq_along(labels))
  colnames(ans)<-c('Group.1',colnames(transposed));ans
}

make_exact_mean_CellChat <- function() {
  fun<-CellChat::computeCommunProb
  body_parts<-as.list(body(fun))
  is_aggregate_assignment<-function(x,name)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name(name))&&is.call(x[[3]])&&identical(x[[3]][[1]],as.name('aggregate'))
  insert<-which(vapply(body_parts,is_aggregate_assignment,logical(1),name='data.use.avg'))
  stopifnot(length(insert)==1L)
  body_parts[[insert]][[3]]<-quote(.cc_exact_group_aggregate(data.use.t,group,FunMean))
  boot<-which(vapply(body_parts,function(x)is.call(x)&&identical(x[[1]],as.name('<-'))&&identical(x[[2]],as.name('data.use.avg.boot')),logical(1)))
  stopifnot(length(boot)==1L)
  call<-body_parts[[boot]][[3]];lambda<-call[['FUN']];nested<-as.list(lambda[[3]])
  j<-which(vapply(nested,is_aggregate_assignment,logical(1),name='data.use.avgB'))
  stopifnot(length(j)==1L)
  nested[[j]][[3]]<-quote(.cc_exact_group_aggregate(data.use.t,groupboot,FunMean))
  lambda[[3]]<-as.call(nested);call[['FUN']]<-lambda;body_parts[[boot]][[3]]<-call
  body(fun)<-as.call(append(body_parts,list(quote(data.use.t<-t(data.use))),after=insert-1L))
  env<-new.env(parent=environment(CellChat::computeCommunProb))
  assign('.cc_exact_group_aggregate',cc_exact_group_aggregate,envir=env);environment(fun)<-env
  fun
}

## Standby only: parallelize the expression means, then keep the much smaller
## per-LR bootstrap calculations in the parent, without changing their formula.
make_exact_mean_serial_LR_CellChat <- function() {
  fun<-make_exact_mean_CellChat();replacements<-0L
  dispatch<-quote(function(X,FUN,...) {
    if(future::nbrOfWorkers()>1 && grepl('.cc_exact_group_aggregate',paste(deparse(body(FUN)),collapse=' '),fixed=TRUE)) {
      future.apply::future_sapply(X=X,FUN=FUN,...)
    } else base::sapply(X=X,FUN=FUN,...)
  })
  replace<-function(expr) {
    if(!is.call(expr))return(expr)
    if(identical(expr[[1]],as.name('<-')) && identical(expr[[2]],as.name('my.sapply'))) {
      replacements<<-replacements+1L;expr[[3]]<-dispatch;return(expr)
    }
    for(i in seq_along(expr))if(i>1 && is.call(expr[[i]]))expr[[i]]<-replace(expr[[i]])
    expr
  }
  body(fun)<-replace(body(fun));stopifnot(replacements==1L)
  fun
}
