## Independent verification of small saved network arrays only, no inference.
args<-commandArgs(TRUE);stopifnot(length(args)%in%c(2L,3L))
allow_pending_roles<-length(args)==3L && identical(args[3],'allow_pending_roles')
out<-args[1];dest<-args[2]
suppressPackageStartupMessages(library(jsonlite))
raw<-readRDS(file.path(out,'network_probabilities_before_significance_filter.rds'))
obj<-readRDS(file.path(out,'CellChat_network.rds'))
prob<-raw$prob;pv<-raw$pval
stopifnot(identical(dim(prob),dim(pv)),all(is.finite(prob)),all(is.finite(pv)),all(prob>=0),all(pv>=0&pv<=1))
filtered<-prob;filtered[pv>=.05]<-0
err<-max(abs(filtered-obj$net$prob));stopifnot(err<1e-12)
count<-apply(filtered>0,c(1,2),sum);weight<-apply(filtered,c(1,2),sum)
cc<-obj$net$count;ww<-obj$net$weight
count_err<-max(abs(cc-count[rownames(cc),colnames(cc),drop=FALSE]))
weight_err<-max(abs(ww-weight[rownames(ww),colnames(ww),drop=FALSE]))
stopifnot(count_err==0,weight_err<1e-12)
db<-obj$DB$interaction;lr<-dimnames(filtered)[[3]]
pw<-db[lr,'pathway_name'];stopifnot(!anyNA(pw))
per_pw_err<-numeric(length(obj$netP$pathways));names(per_pw_err)<-obj$netP$pathways
for(pathway in names(per_pw_err)) {
 expected<-apply(filtered[,,pw==pathway,drop=FALSE],c(1,2),sum)
 actual<-obj$netP$prob[,,pathway,drop=FALSE][,,1]
 per_pw_err[pathway]<-max(abs(expected-actual))
}
stopifnot(all(per_pw_err<1e-12))
tab<-read.delim(file.path(out,'ligand_receptor_communications.tsv'),check.names=FALSE)
stopifnot(nrow(tab)==sum(filtered>0),all(tab$pval<.05),all(tab$prob>0))
role_path<-file.path(out,'cellchat_role_scatter_data.tsv')
role_err<-NULL
if(file.exists(role_path)) {
 role<-read.delim(role_path,check.names=FALSE)
 role_err<-numeric(nrow(role))
 for(i in seq_len(nrow(role))) {
  q<-role$signaling[i];ccell<-role$cell_type[i]
  if(q=='All')mat<-count else mat<-apply(filtered[,,pw%in%q,drop=FALSE]>0,c(1,2),sum)
  expected<-rowSums(mat)+colSums(mat)-diag(mat)
  role_err[i]<-abs(role$LR_count[i]-expected[ccell])
 }
 stopifnot(all(is.finite(role_err)),all(role_err==0))
} else {
 stopifnot(allow_pending_roles)
}
groups<-read.delim(file.path(out,'modeled_group_counts.tsv'),check.names=FALSE)
stopifnot(nrow(groups)==8,sum(groups$cell_count)==528739)
report<-list(status=if(is.null(role_err))'preliminary_saved_arrays_pending_rendering' else 'validated_saved_network_scope',scope='Only small saved probabilities, pathway sums, network tables and role counts when available; no algorithms rerun',
 modeled_records=sum(groups$cell_count),modeled_groups=nrow(groups),RNA_records_retained=581725,
 pvalue_cutoff='Strict nominal p < .05; no donor-level or FDR claim',
 probability_max_error_after_strict_filter=err,network_count_max_error=count_err,network_weight_max_error=weight_err,
 pathway_sum_max_error=max(per_pw_err),role_scatter_LR_count_max_error=if(is.null(role_err))NULL else max(role_err),
 total_LR_records=nrow(tab),unique_LR=length(unique(tab$interaction_name)),pathways=length(per_pw_err),
 prefilter_positive_prob_p_equals_005=sum(prob>0 & pv==.05),all_p_equals_005_edges_excluded=all(obj$net$prob[pv==.05]==0),
 count_matrix_groups=rownames(cc),statistical_scope='Transcript-based pooled mixed-cohort hypotheses; neither actual measured contact nor AD-specific or AD-HC comparison')
write_json(report,dest,pretty=TRUE,auto_unbox=TRUE,digits=17)
cat(report$status,'\n')
