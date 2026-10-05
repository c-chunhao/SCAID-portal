# Process-local interface adapter for scDblFinder 1.16.0 + xgboost 3.2.1.1.
# Preserve the original scDblFinder classifier design while using current
# DMatrix/params/evals contracts. Never modify an installed package on disk.
.adapter_audit <- new.env(parent=emptyenv())
.adapter_audit$calls <- list(); .adapter_audit$errors <- list()
.adapter_audit$channel <- '';.adapter_audit$probe <- NULL;.adapter_audit$previous <- NULL
.ns_scDbl <- asNamespace('scDblFinder')
.original_xgbtrain <- get('.xgbtrain',.ns_scDbl)
.modern_xgbtrain <- function(d2,ctype,nrounds=NULL,max_depth=6,nfold=5,tree_method='exact',subsample=.75,nthreads=1,metric='logloss',...) {
 call_id <- length(.adapter_audit$calls)+1L
 tryCatch({
  if(!is.integer(ctype))ctype<-as.integer(ctype)-1L
  d2<-as.matrix(d2);stopifnot(length(ctype)==nrow(d2),all(ctype %in% c(0L,1L)),length(unique(ctype))==2L)
  dtrain<-xgboost::xgb.DMatrix(data=d2,label=ctype)
  stopifnot(identical(as.numeric(xgboost::getinfo(dtrain,'label')),as.numeric(ctype)))
  extra<-list(...)
  train_params<-modifyList(list(max_depth=max_depth,objective='binary:logistic',eval_metric=metric,tree_method=tree_method,nthread=nthreads),extra)
  cv_params<-modifyList(train_params,list(subsample=subsample))
  if(is.null(nrounds))nrounds<-0L
  stopifnot(is.numeric(nrounds),nrounds>=0)
  cv_ran<-FALSE;cv_best<-NA_integer_
  if(nrounds<=1) {
   cv_ran<-TRUE
   cv<-xgboost::xgb.cv(params=cv_params,data=dtrain,nrounds=200,nfold=nfold,early_stopping_rounds=2,verbose=FALSE)
   stopifnot(!is.null(cv$evaluation_log),nrow(cv$evaluation_log)>0)
   e<-cv$evaluation_log
   mean_name<-grep('test.+mean',colnames(e),value=TRUE)[1]
   std_name<-grep('test.+std',colnames(e),value=TRUE)[1]
   # The official modern return records early stopping via callbacks. Choosing
   # its observed minimum test metric also makes the selected round explicit.
   cv_best<-if(!is.null(cv$best_iteration))as.integer(cv$best_iteration) else which.min(e[[mean_name]])
   if(nrounds==0)nrounds<-cv_best else {
    bound<-e[[mean_name]][cv_best]+nrounds*e[[std_name]][cv_best]
    nrounds<-min(which(e[[mean_name]]<=bound))
   }
  }
  stopifnot(nrounds>=1)
  # Legacy high-level xgboost used the training set as its evaluation watchlist.
  fit<-xgboost::xgb.train(params=train_params,data=dtrain,nrounds=as.integer(nrounds),evals=list(train=dtrain),early_stopping_rounds=2,verbose=FALSE)
  stopifnot(inherits(fit,'xgb.Booster'))
  trained_rounds<-xgboost::xgb.get.num.boosted.rounds(fit)
  stopifnot(trained_rounds>=1L,trained_rounds<=nrounds)
  if(is.null(.adapter_audit$probe)).adapter_audit$probe<-d2
  prob<-predict(fit,.adapter_audit$probe)
  stopifnot(is.numeric(prob),all(is.finite(prob)),all(prob>=0&prob<=1))
  delta<-if(is.null(.adapter_audit$previous))NA_real_ else mean(abs(prob-.adapter_audit$previous))
  .adapter_audit$previous<-prob
  .adapter_audit$calls[[call_id]]<-data.frame(call_id=call_id,channel=.adapter_audit$channel,training_rows=nrow(d2),positive_labels=sum(ctype==1L),negative_labels=sum(ctype==0L),DMatrix_labels_exact=TRUE,CV_executed=cv_ran,CV_best_round=cv_best,train_rounds=as.integer(nrounds),actual_trained_rounds=trained_rounds,probe_score_sd=sd(prob),probe_mean_change_from_previous=delta,probe_score_SHA256=digest::digest(prob,algo='sha256'),stringsAsFactors=FALSE)
  fit
 },error=function(e){
  .adapter_audit$errors[[length(.adapter_audit$errors)+1L]]<-paste(.adapter_audit$channel,conditionMessage(e))
  stop(e)
 })
}
unlockBinding('.xgbtrain',.ns_scDbl)
assign('.xgbtrain',.modern_xgbtrain,envir=.ns_scDbl)
lockBinding('.xgbtrain',.ns_scDbl)

# scDblFinder 1.16.0 also catches errors from predict(), outside .xgbtrain.
# Replace just that handler in this process so every classifier failure is
# explicit, including a failure after a successful Booster fit.
.original_scDblscore <- get('.scDblscore',.ns_scDbl)
.strict_scDblscore <- .original_scDblscore
.strict_handler_count <- 0L
.no_fallback_handler <- quote(function(e) {
 audit <- get('.adapter_audit',envir=.GlobalEnv)
 audit$errors[[length(audit$errors)+1L]] <- paste(audit$channel,'upstream classifier:',conditionMessage(e))
 stop(e)
})
.rewrite_classifier_error_handler <- function(expr) {
 if(!is.call(expr))return(expr)
 expr <- as.call(lapply(as.list(expr),.rewrite_classifier_error_handler))
 if(identical(expr[[1]],as.name('tryCatch')) && 'error' %in% names(expr)) {
  old <- paste(deparse(expr[['error']]),collapse=' ')
  if(grepl('function.*\\(e\\).*d\\$score',old)) {
   expr[['error']] <- .no_fallback_handler
   .strict_handler_count <<- .strict_handler_count+1L
  }
 }
 expr
}
body(.strict_scDblscore) <- .rewrite_classifier_error_handler(body(.strict_scDblscore))
stopifnot(.strict_handler_count==1L)
unlockBinding('.scDblscore',.ns_scDbl)
assign('.scDblscore',.strict_scDblscore,envir=.ns_scDbl)
lockBinding('.scDblscore',.ns_scDbl)
