## Read original objects sequentially; never mutate original RDS or live database.
suppressPackageStartupMessages({library(CellChat); library(Matrix); library(digest)})
root <- Sys.getenv('SCAID_CELLCHAT_AUDIT_DIR', unset='')
if (!nzchar(root)) stop('Set SCAID_CELLCHAT_AUDIT_DIR to the local audit output/input directory')
future::plan('sequential'); set.seed(20261005)
inv <- read.delim(file.path(root, 'tables/cellchat_explicit_inventory.tsv'), check.names=FALSE)
subset_pattern <- Sys.getenv('CELLCHAT_AUDIT_SUBSET','')
if(nzchar(subset_pattern)) inv <- inv[grepl(subset_pattern,inv$cellchat_path),]
suffix <- Sys.getenv('CELLCHAT_AUDIT_SUFFIX','')
outfile <- function(name)file.path(root,'tables',sub('(\\.[^.]+)$',paste0(suffix,'\\1'),name))
summaries <- list(); groups <- list(); pathways <- list(); parameters <- list()
for (i in seq_len(nrow(inv))) {
  f <- inv$cellchat_path[i]; id <- if('artifact_id' %in% names(inv))inv$artifact_id[i]else inv$object_stem[i]
  message(sprintf('[%d/%d] %s', i, nrow(inv), id))
  tryCatch({
    cc <- readRDS(f)
    stopifnot(inherits(cc,'CellChat'))
    meta <- cc@meta; ident <- as.character(cc@idents)
    pr <- cc@net$prob; pv <- cc@net$pval
    stopifnot(length(ident)==nrow(meta), identical(dim(pr),dim(pv)), length(dim(pr))==3L)
    groupnames <- dimnames(pr)[[1]]
    stopifnot(identical(groupnames,dimnames(pr)[[2]]),setequal(groupnames,unique(ident)))
    n <- table(factor(ident,levels=groupnames))
    sample_col <- intersect(c('Sample','sample','orig.ident'),names(meta))[1]
    disease_col <- intersect(c('Disease','disease','condition'),names(meta))[1]
    donor_col <- intersect(c('ind_cov','donor_id','donor','individual'),names(meta))[1]
    ng <- length(groupnames); nl <- dim(pr)[3]
    finite <- all(is.finite(pr)) && all(is.finite(pv)); nonnegative <- all(pr>=0)
    valid <- pv < 0.05 & pr > 0
    dbhash <- digest(cc@DB,algo='sha256',serialize=TRUE)
    dbn <- nrow(cc@DB$interaction)
    lrnames <- dimnames(pr)[[3]]
    map_idx <- match(lrnames, rownames(cc@DB$interaction))
    if (all(is.na(map_idx)) && 'interaction_name' %in% names(cc@DB$interaction)) map_idx <- match(lrnames,cc@DB$interaction$interaction_name)
    lrpath <- as.character(cc@DB$interaction$pathway_name[map_idx])
    # Original permutation p-values are conditional on this pooled object.
    # Export no disease p-values; raw probabilities are not observed communication.
    for (path in sort(unique(lrpath[!is.na(lrpath)]))) {
      ix <- which(lrpath==path)
      prob <- pr[,,ix,drop=FALSE]; sig <- valid[,,ix,drop=FALSE]
      # Distinguish evaluated zero from unknown/unmatched pathways.
      pathways[[length(pathways)+1L]] <- data.frame(object_stem=id,pathway=path,evaluated_ligand_receptor_pairs=length(ix),cell_groups=ng,potential_group_pairs=ng^2,significant_edges=sum(sig),sum_conditional_probability=sum(prob[sig]),mean_probability_per_evaluated_group_LR_pair=sum(prob[sig])/(ng^2*length(ix)),scope='pooled-object descriptive; not a donor disease effect')
    }
    for (g in groupnames) {
      gix <- ident==g
      groups[[length(groups)+1L]] <- data.frame(object_stem=id,cell_group=g,n_cells=unname(n[g]),n_specimen_labels=if(!is.na(sample_col))length(unique(meta[[sample_col]][gix]))else NA_integer_)
    }
    summaries[[length(summaries)+1L]] <- data.frame(object_stem=id,path=f,n_cells=nrow(meta),n_groups=ng,min_group_cells=min(n),n_specimen_labels=if(!is.na(sample_col))length(unique(meta[[sample_col]]))else NA_integer_,n_verified_individual_keys=if(!is.na(donor_col))length(unique(meta[[donor_col]]))else NA_integer_,disease_values=if(!is.na(disease_col))paste(sort(unique(as.character(meta[[disease_col]]))),collapse=';')else 'not_saved',sample_column=if(!is.na(sample_col))sample_col else '',donor_column=if(!is.na(donor_col))donor_col else '',n_lr_tested=nl,lr_database_rows=dbn,lr_database_SHA256=dbhash,n_lr_missing_DB_map=sum(is.na(map_idx)),probability_finite=finite,probability_nonnegative=nonnegative,pvalues_in_range=all(pv>=0 & pv<=1),network_dimensions_match=TRUE,n_p05_positive_edges=sum(valid),all_group_comparisons=ng^2*nl,reported_package_version=as.character(packageVersion('CellChat')),source_package_version='not recorded in old object',minimum_group_lt10=sum(n<10),status='audited original; descriptive only')
    parameters[[id]] <- cc@options
    saveRDS(cc@DB,file.path(root,'tables',paste0('CellChatDB_',dbhash,'.rds')))
    write.table(do.call(rbind,summaries),outfile('cellchat_object_QA.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    write.table(do.call(rbind,groups),outfile('cellchat_group_cell_counts.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    write.table(do.call(rbind,pathways),outfile('cellchat_evaluated_pathway_descriptive.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
    saveRDS(parameters,outfile('cellchat_saved_options.rds'))
    rm(cc,meta,pr,pv,valid);gc()
  },error=function(e){
    write(paste(id,conditionMessage(e),sep='\t'),file.path(root,'logs/cellchat_read_errors.tsv'),append=TRUE)
    gc()
  })
}
writeLines(capture.output(sessionInfo()),file.path(root,'logs/cellchat_R_sessionInfo.txt'))
message('Completed explicit CellChat object audit.')
