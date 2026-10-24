## Independent consistency check of the saved network, tables and cutoff.
args<-commandArgs(TRUE);stopifnot(length(args)==2L)
out<-args[1];formal<-args[2]=='formal'
suppressPackageStartupMessages({library(jsonlite);library(CellChat)})
network<-readRDS(file.path(out,'CellChat_network.rds'))
raw<-readRDS(file.path(out,'network_probabilities_before_significance_filter.rds'))
computed<-read_json(file.path(out,'cellchat_computation.json'),simplifyVector=TRUE)
stopifnot(identical(dim(raw$prob),dim(network$net$prob)),identical(dim(raw$pval),dim(network$net$pval)))
expected<-raw$prob;expected[raw$pval>=.05]<-0
stopifnot(identical(expected,network$net$prob),all(network$net$prob[raw$pval>=.05]==0))
expected_count<-apply(expected>0,c(1,2),sum)
expected_weight<-apply(expected,c(1,2),sum)
stopifnot(max(abs(expected_count-network$net$count))==0,max(abs(expected_weight-network$net$weight))<1e-12)
for(which in c('count','weight')) {
  tab<-read.delim(file.path(out,paste0('network_',which,'.tsv')),check.names=FALSE,stringsAsFactors=FALSE)
  values<-as.matrix(tab[,-1,drop=FALSE]);rownames(values)<-tab$source
  stopifnot(identical(rownames(values),rownames(network$net[[which]])),identical(colnames(values),colnames(network$net[[which]])),max(abs(values-network$net[[which]]))<1e-12)
}
lr<-read.delim(file.path(out,'ligand_receptor_communications.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(lr$prob>0),all(lr$pval<.05),nrow(lr)==sum(expected_count),nrow(lr)==computed$significant_communication_records)
groups<-read.delim(file.path(out,'modeled_group_counts.tsv'),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(sum(groups$cell_count)==computed$modeled_cells,nrow(groups)==length(computed$modeled_cell_types),all(groups$cell_count>=10))
if(formal) {
  stopifnot(computed$input_cells==581725,computed$modeled_cells==528739,nrow(groups)==8,
    computed$excluded_ambiguous==52809,computed$excluded_provisional_nonimmune==172,
    computed$excluded_small_group==5,computed$nboot==100,computed$DB_genes_in_RNA==1434,
    computed$restored_RNA_symbol_rows==9)
}
answer<-list(status='validated',formal=formal,input_cells=computed$input_cells,modeled_cells=computed$modeled_cells,
  groups=nrow(groups),bootstrap=computed$nboot,strict_p_less_than_005=TRUE,
  saved_probabilities_match_unfiltered_probability_mask_exactly=TRUE,
  network_count_max_absolute_error=max(abs(expected_count-network$net$count)),
  network_weight_max_absolute_error=max(abs(expected_weight-network$net$weight)),
  significant_LR_records=nrow(lr),unique_LR_interactions=length(unique(lr$interaction_name)),
  all_network_tables_match_RDS=TRUE,source_scope=computed$scope)
write_json(answer,file.path(out,'independent_network_consistency.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat('FINAL_CELLCHAT_NETWORK_CONSISTENCY_VALIDATED\n')
