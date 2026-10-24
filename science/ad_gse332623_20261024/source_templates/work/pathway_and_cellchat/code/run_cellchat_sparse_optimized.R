## Standby internal implementation optimization; original R script unchanged.
args<-commandArgs(TRUE);cfg<-jsonlite::fromJSON(args[1])
source(file.path(cfg$code_dir,'pathway_core.R'))
source(file.path(cfg$code_dir,'symbol_alias_view.R'))
source(file.path(cfg$code_dir,'cellchat_exact_mean.R'))
source(cfg$sparse_helper_R)
load_sparse_trimean_prototype(cfg$sparse_helper_cpp,file.path(cfg$output,'compiled_cache'))
suppressPackageStartupMessages({library(CellChat);library(future);library(ggplot2);library(ComplexHeatmap)})
if(cfg$workers>1)future::plan('multisession',workers=cfg$workers) else future::plan('sequential')
options(future.globals.maxSize=80*1024^3)
worker_pids<-if(cfg$workers>1)unique(unlist(future.apply::future_lapply(seq_len(cfg$workers),function(i)Sys.getpid(),future.scheduling=1))) else integer()
jsonlite::write_json(list(parent_R_pid=Sys.getpid(),worker_pids=worker_pids),file.path(cfg$output,'cellchat_processes.json'),auto_unbox=FALSE,pretty=TRUE)
set.seed(cfg$seed)
out<-cfg$output;dir.create(out,recursive=TRUE,showWarnings=FALSE)
axes<-read_axes(cfg$input)
meta<-read.delim(file.path(out,'cellsmeta.tsv'),check.names=FALSE,stringsAsFactors=FALSE,quote='"')
stopifnot(identical(meta$cell_id,axes$cells))
labels<-trimws(as.character(meta$MajorCellType))
unknown<-is.na(labels)|!nzchar(labels)|grepl('^(unknown|unknonw|unassigned|unresolved|ambiguous|mixed|doublet|other|others)(\\b|[ _/:-])',labels,ignore.case=TRUE)
provisional<-rep(FALSE,length(labels))
if(!is.null(cfg$coarse_review_path)) {
  review<-read.delim(cfg$coarse_review_path,check.names=FALSE,stringsAsFactors=FALSE)
  provisional_types<-review$MajorCellType[grepl('source_dominated',review$review_status,ignore.case=TRUE)]
  provisional<-labels%in%provisional_types
}
counts<-table(labels[!unknown & !provisional]);small<-names(counts)[counts<cfg$min_cells]
keep<-!unknown & !provisional & !labels%in%small
reason<-ifelse(unknown,'Unresolved or ambiguous annotation',ifelse(provisional,'source_dominated_provisional_nonimmune_prediction',ifelse(labels%in%small,paste0('Fewer than ',cfg$min_cells,' cells'),'Included')))
included<-data.frame(cell_id=axes$cells,MajorCellType=labels,included=keep,reason=reason)
write.table(included,gzfile(file.path(out,'cell_inclusion.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE)
stopifnot(sum(keep)>0,length(unique(labels[keep]))>=2)
db<-CellChatDB.human;saveRDS(db,file.path(out,'CellChatDB_human_used.rds'))
db_genes<-unique(extractGene(db))
view<-if(!is.null(cfg$symbol_alias_path))canonical_rna_view(axes$genes,cfg$symbol_alias_path) else list(genes=axes$genes)
gene_idx<-which(view$genes%in%db_genes);genes<-view$genes[gene_idx];source_genes<-axes$genes[gene_idx]
write.table(data.frame(DB_gene=db_genes,present_in_RNA=db_genes%in%genes,
  source_RNA_symbol=ifelse(db_genes%in%genes,source_genes[match(db_genes,genes)],NA_character_)),file.path(out,'CellChatDB_gene_overlap.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(canonical_symbol=genes,source_symbol=source_genes,alias_restored=genes!=source_genes),file.path(out,'CellChat_RNA_symbol_matching_view.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
## Internal performance cache: exact DB-gene columns from the already
## verified source-RNA CSC index, before the same cell-inclusion mask.
cache<-H5File$new(cfg$DB_expression_cache,'r')
cache_cells<-as.character(cache[['obs/cell_id']][])
cache_sources<-as.character(cache[['gene/source_symbols']][])
cache_canonical<-as.character(cache[['gene/canonical_symbols']][])
stopifnot(identical(cache_cells,axes$cells),identical(cache_sources,source_genes),identical(cache_canonical,genes))
cache_keep<-as.logical(cache[['obs/include_cellchat']][])
stopifnot(identical(cache_keep,keep))
x<-cache[['expression_csc/data']][]
i<-as.integer(cache[['expression_csc/indices']][])
p<-cache[['expression_csc/indptr']][]
stopifnot(length(x)==length(i),length(p)==length(genes)+1L,max(p)<.Machine$integer.max)
mc<-new('dgCMatrix',i=i,x=as.numeric(x),p=as.integer(p),Dim=as.integer(c(length(axes$cells),length(genes))),Dimnames=list(axes$cells,genes))
cache$close_all();rm(x,i,p,cache_cells,cache_sources,cache_canonical,cache_keep);gc(FALSE)
m<-as(t(mc),'dgCMatrix')[,keep,drop=FALSE]
rm(mc);gc(FALSE)
stopifnot(identical(colnames(m),axes$cells[keep]),all(is.finite(m@x)),all(m@x>=0))
ccmeta<-data.frame(MajorCellType=factor(labels[keep],levels=sort(unique(labels[keep]))),row.names=axes$cells[keep])
ccmeta$samples<-factor(paste0(gsub('/','_',cfg$object_id),'_mixed_cohort'))
clock<-proc.time()[['elapsed']]
cat(jsonlite::toJSON(list(stage='matrix_ready',RNA_genes=length(axes$genes),DB_genes=length(genes),cells=ncol(m),groups=nlevels(ccmeta$MajorCellType),RSS_GiB=current_rss_gib()),auto_unbox=TRUE),'\n');flush.console()
cellchat<-createCellChat(object=m,meta=ccmeta,group.by='MajorCellType',do.sparse=TRUE)
rm(m);gc(FALSE);cellchat@DB<-db
cellchat<-subsetData(cellchat)
cellchat<-identifyOverExpressedGenes(cellchat,do.fast=TRUE,min.cells=cfg$min_cells,only.pos=TRUE,thresh.pc=0,thresh.fc=0,thresh.p=0.05)
cellchat<-identifyOverExpressedInteractions(cellchat)
equivalent_computeCommunProb<-make_positive_order_CellChat()
## Progress logging wraps the same pure helper without changing its inputs,
## return value, group order or random-number state.
helper_env<-environment(equivalent_computeCommunProb)
progress_helper<-local({
  original_helper<-get('.cc_positive_group_aggregate',envir=helper_env)
  completed<-0L
  function(transposed,group,FUN) {
    value<-original_helper(transposed,group,FUN)
    completed<<-completed+1L
    cat(jsonlite::toJSON(list(stage='sparse_exact_triMean_completed',completed=completed,total=cfg$nboot+1L,utc=as.character(Sys.time())),auto_unbox=TRUE),'\n')
    flush.console();value
  }
})
assign('.cc_positive_group_aggregate',progress_helper,envir=helper_env)
cellchat<-equivalent_computeCommunProb(cellchat,type='triMean',raw.use=TRUE,population.size=FALSE,distance.use=FALSE,nboot=cfg$nboot,seed.use=cfg$seed)
cellchat<-filterCommunication(cellchat,min.cells=cfg$min_cells)
saveRDS(list(prob=cellchat@net$prob,pval=cellchat@net$pval),file.path(out,'network_probabilities_before_significance_filter.rds'),compress='gzip')
cellchat@net$prob[cellchat@net$pval>=0.05]<-0
cellchat<-computeCommunProbPathway(cellchat,thresh=0.05)
cellchat<-aggregateNet(cellchat,thresh=0.05)
## Some igraph centrality routines initialize randomly; keep that small step
## sequential and explicitly seeded rather than using unseeded future workers.
future::plan('sequential');set.seed(cfg$seed)
cellchat<-netAnalysis_computeCentrality(cellchat,slot.name='netP')
saveRDS(cellchat,file.path(out,'CellChat_corrected.rds'),compress='gzip')
saveRDS(list(net=cellchat@net,netP=cellchat@netP,idents=cellchat@idents,DB=cellchat@DB,options=cellchat@options),file.path(out,'CellChat_network.rds'),compress='gzip')
df<-subsetCommunication(cellchat)
write.table(df,file.path(out,'ligand_receptor_communications.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
write.table(subsetCommunication(cellchat,slot.name='netP'),file.path(out,'pathway_communications.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
for(name in c('count','weight'))write.table(data.frame(source=rownames(cellchat@net[[name]]),cellchat@net[[name]],check.names=FALSE),file.path(out,paste0('network_',name,'.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
write.table(data.frame(MajorCellType=names(table(cellchat@idents)),cell_count=as.integer(table(cellchat@idents))),file.path(out,'modeled_group_counts.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
report<-list(status='computed',object_id=cfg$object_id,mode=cfg$mode,input_cells=length(axes$cells),modeled_cells=sum(keep),excluded_ambiguous=sum(unknown),excluded_provisional_nonimmune=sum(provisional & !unknown),excluded_small_group=sum(!unknown & !provisional & !keep),min_cells=cfg$min_cells,
  modeled_cell_types=levels(cellchat@idents),human_DB_interactions=nrow(db$interaction),RNA_genes=length(axes$genes),DB_genes_in_RNA=length(genes),
  significant_communication_records=nrow(df),unique_LR_interactions=length(unique(df$interaction_name)),inferred_pathways=cellchat@netP$pathways,nboot=cfg$nboot,seed=cfg$seed,workers=cfg$workers,centrality_workers=1,centrality_seed=cfg$seed,type='triMean',population_size=FALSE,raw_use=TRUE,distance_use=FALSE,edge_significance='p < 0.05; unfiltered probabilities retained separately',
  expression_input='RNA log1p(cp10k), originally normalized on all RNA genes; DB subset not renormalized',PPI_projection='Not applied; current CellChat standard raw.use=TRUE workflow',
  symbol_alias_mapping_sha256=cfg$symbol_alias_sha256,symbol_alias_mapping=cfg$symbol_alias_path,
  restored_RNA_symbol_rows=sum(genes!=source_genes),symbol_alias_scope='Nine verified one-to-one source Ensembl/local annotation mappings for internal DB matching; original RNA axis and values unchanged',
  implementation='Exact nonnegative sparse-positive quartiles using the installed collapse type-7 arithmetic; official R mean of Q25,Q50,Q50,Q75 retained. Single-worker means and sequential LR bootstrap; full 100-bootstrap pilot and full-axis aggregation equivalence independently checked.',
  DB_expression_cache_sha256=cfg$DB_expression_cache_sha256,DB_expression_cache_receipt_sha256=cfg$DB_expression_cache_receipt_sha256,implementation_equivalence_sha256=cfg$implementation_equivalence_sha256,
  implementation_review_sha256=cfg$implementation_review_sha256,full_axis_aggregation_review_sha256=cfg$full_axis_aggregation_review_sha256,installed_official_function_source_sha256=cfg$installed_official_function_source_sha256,
  scope='Predicted communication among annotated cell populations pooled across the mixed study cohort. No AD-specific or AD-HC comparison is inferred.',
  role_scatter_workaround='Compute point Count directly from significant LR probabilities; avoid aggregateNet signaling-filter source-target regex splitting bug in CellChat 2.2.0.9001',
  seconds=proc.time()[['elapsed']]-clock,RSS_GiB=current_rss_gib(),packages=as.list(sapply(c('CellChat','Matrix','hdf5r','future','presto'),function(x)as.character(packageVersion(x)))))
write_json(report,file.path(out,'cellchat_computation.json'),auto_unbox=TRUE,pretty=TRUE,digits=17)
cat('CELLCHAT_COMPUTED\n')
source(file.path(cfg$code_dir,'render_cellchat.R'))
render_cellchat(cellchat,out,cfg)
