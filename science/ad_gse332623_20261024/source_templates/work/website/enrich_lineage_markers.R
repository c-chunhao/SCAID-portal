## Recomputed lineage-marker enrichment, not AD-versus-HC enrichment.
args <- commandArgs(TRUE)
stopifnot(length(args) %in% c(3L,4L))
plots_only <- length(args)==4L && args[4]=="--plots-only"
stopifnot(length(args)==3L || plots_only)
rna_h5 <- args[1]; marker_tsv <- args[2]; out <- args[3]
.libPaths(c('@R_LIB@', .libPaths()))
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db); library(rhdf5)
  library(ggplot2); library(jsonlite)
})
dir.create(out, recursive=TRUE, showWarnings=FALSE)
sha <- function(p) substr(system2('sha256sum', shQuote(p), stdout=TRUE), 1, 64)
if (!plots_only) {
var_index <- h5readAttributes(rna_h5, 'var')[['_index']]
genes <- as.character(h5read(rna_h5, paste0('var/', var_index)))
stopifnot(!anyNA(genes), !anyDuplicated(genes))
source_genes <- genes
alias_file <- Sys.getenv('SCAID_VERIFIED_ALIAS_TABLE',
  file.path(dirname(dirname(rna_h5)),'provenance','resource_gene_aliases.tsv'))
stopifnot(file.exists(alias_file))
aliases <- read.delim(alias_file,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(all(c('source_symbol','canonical_symbol','source_ensembl_id') %in% names(aliases)),
          !anyDuplicated(aliases$source_symbol), !anyDuplicated(aliases$canonical_symbol),
          all(aliases$source_symbol %in% source_genes),
          !any(aliases$canonical_symbol %in% source_genes))
canonicalize <- function(x) {
  location <- match(x,aliases$source_symbol); mapped <- !is.na(location)
  x[mapped] <- aliases$canonical_symbol[location[mapped]];x
}
genes <- canonicalize(source_genes)
stopifnot(!anyDuplicated(genes))
write.table(aliases,file.path(out,'resource_gene_aliases.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
listing <- h5ls(rna_h5, recursive=TRUE)
index_row <- listing[listing$group=='/X' & listing$name=='indices',,drop=FALSE]
stopifnot(nrow(index_row)==1L)
nnz <- as.numeric(strsplit(index_row$dim, ' x ', fixed=TRUE)[[1]][1])
detected <- numeric(length(genes))
for (start in seq(1, nnz, by=2000000)) {
  values <- h5read(rna_h5, 'X/indices', start=start, count=min(2000000, nnz-start+1))
  detected <- detected + tabulate(as.integer(values)+1L, nbins=length(genes))
}
universe <- genes[detected>=3]
write.table(data.frame(gene=source_genes, annotation_symbol=genes, detected_cells=detected, enrichment_background=detected>=3),
            file.path(out,'gene_background.tsv'), sep='\t', quote=FALSE, row.names=FALSE)
mk <- read.delim(marker_tsv, check.names=FALSE, stringsAsFactors=FALSE)
mk$source_gene <- mk$gene
mk$gene <- canonicalize(mk$gene)
required <- c('p_val_adj','avg_log2FC','pct.1','gene','group1','group2')
stopifnot(all(required %in% names(mk)), all(mk$group2=='others'))
good <- is.finite(mk$p_val_adj) & mk$p_val_adj<.05 &
        is.finite(mk$avg_log2FC) & mk$avg_log2FC>1 & mk[['pct.1']]>=.1
mk <- mk[good & mk$gene %in% universe & !mk$group1 %in% c('Unassigned','Unknown'),,drop=FALSE]
sets_file <- '@PROJECT@/20240526_kegg_hsa_list.rds'
sets <- readRDS(sets_file)
t2g <- do.call(rbind, lapply(names(sets), function(id)
  data.frame(term=substr(id,1,8), gene=unique(sets[[id]]), stringsAsFactors=FALSE)))
t2g$gene <- canonicalize(t2g$gene)
t2g <- unique(t2g)
t2n <- data.frame(term=substr(names(sets),1,8),
                  name=gsub('_',' ',substring(names(sets),10)), stringsAsFactors=FALSE)
results <- list(GO=list(), KEGG=list()); eligibility <- list()
for (group in unique(mk$group1)) {
  selected <- unique(mk$gene[mk$group1==group])
  eligibility[[group]] <- length(selected)
  if (length(selected)<5L) next
  message('ENRICH ', group, ' markers=', length(selected))
  go <- enrichGO(selected, universe=universe, OrgDb=org.Hs.eg.db, keyType='SYMBOL',
                 ont='BP', pvalueCutoff=1, qvalueCutoff=1, pAdjustMethod='BH',
                 minGSSize=10, maxGSSize=500, readable=FALSE)
  kg <- enricher(selected, universe=universe, TERM2GENE=t2g, TERM2NAME=t2n,
                 pvalueCutoff=1, qvalueCutoff=1, pAdjustMethod='BH',
                 minGSSize=10, maxGSSize=500)
  for (method in c('GO','KEGG')) {
    frame <- as.data.frame(if(method=='GO') go else kg)
    if (nrow(frame)) { frame$cell_type <- group; results[[method]][[group]] <- frame }
  }
}
} else {
  previous <- read_json(file.path(out,'enrichment_receipt.json'),simplifyVector=TRUE)
  ready <- read_json(file.path(dirname(rna_h5),'analysis_ready.json'),simplifyVector=TRUE)
  stopifnot(previous$input_RNA_sha256==ready$RNA_sha256,
            previous$marker_sha256==sha(marker_tsv))
  results <- lapply(c('GO','KEGG'),function(method){
    frame <- read.delim(file.path(out,paste0(method,'_all_tested.tsv')),check.names=FALSE)
    split(frame,frame$cell_type)
  }); names(results)<-c('GO','KEGG')
}
plot_counts <- list(); effective_backgrounds <- list()
for (method in names(results)) {
  frame <- if (length(results[[method]])) do.call(rbind, results[[method]]) else data.frame()
  if (!plots_only) write.table(frame, file.path(out,paste0(method,'_all_tested.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  bg <- unique(as.integer(sub('^.*/','',frame$BgRatio)))
  stopifnot(length(bg)==1L, !is.na(bg)); effective_backgrounds[[method]]<-bg
  sig <- if(nrow(frame))frame[is.finite(frame$p.adjust)&frame$p.adjust<.05,,drop=FALSE]else frame
  if (!nrow(sig)) { plot_counts[[method]] <- 0L; next }
  top <- do.call(rbind,lapply(split(sig,sig$cell_type),function(d)head(d[order(d$p.adjust,-d$Count),],3)))
  ids <- unique(top$ID); show <- sig[sig$ID %in% ids,,drop=FALSE]
  labels <- setNames(vapply(unique(show$ID),function(id){
    z<-show$Description[match(id,show$ID)];paste(strwrap(z,width=55),collapse='\n')},character(1)),unique(show$ID))
  show$term <- factor(show$ID,levels=rev(ids),labels=labels[rev(ids)])
  show$strength <- pmin(-log10(pmax(show$p.adjust,1e-300)),20)
  title <- if(method=='GO')'GO biological processes associated with cell-type markers'else'KEGG pathways associated with cell-type markers'
  p <- ggplot(show,aes(x=cell_type,y=term,size=Count,color=strength))+
    geom_point(alpha=.9)+scale_color_gradient(low='#8ecae6',high='#a71930',name='-log10 adjusted P\n(capped at 20)')+
    scale_size_continuous(range=c(2,7),name='Marker genes')+
    scale_x_discrete(labels=function(x)ifelse(x=='Fibroblasts','Fibroblasts (provisional)',x))+
    labs(x=NULL,y=NULL,title=title)+theme_bw(base_size=10)+
    theme(panel.grid.minor=element_blank(),axis.text.x=element_text(angle=45,hjust=1),plot.title=element_text(size=11))
  height <- max(5, length(ids)*.40+1.6)
  ggsave(file.path(out,paste0('MajorCellType_',if(method=='GO')'Go'else'KEGG','_Ht.png')),p,width=13,height=height,dpi=160,limitsize=FALSE)
  ggsave(file.path(out,paste0('MajorCellType_',if(method=='GO')'Go'else'KEGG','_Ht.pdf')),p,width=13,height=height,limitsize=FALSE)
  write.table(show,file.path(out,paste0(method,'_displayed.tsv')),sep='\t',quote=FALSE,row.names=FALSE)
  plot_counts[[method]] <- length(ids)
}
record <- if (plots_only) previous else list(status='completed',comparison='cell-type markers within mixed cohort; no AD/HC test',
                input_RNA_sha256=sha(rna_h5),marker_sha256=sha(marker_tsv),
                gene_set_sha256=sha(sets_file),background_genes=length(universe),
                resource_gene_alias_table_sha256=sha(alias_file),
                resource_gene_alias_count=nrow(aliases),
                gene_identifier_rule='Ensembl-backed, collision-free source-symbol to canonical resource-key matching only; source RNA axis and expression values unchanged',
                background_rule='RNA symbol detected in at least 3 QC-retained cells',
                marker_rule='positive log2FC>1, BH adjusted P<0.05, detection>=0.1; Unassigned excluded',
                selected_marker_counts=eligibility,GO_ontology='BP',
                min_gene_set_size=10,max_gene_set_size=500,displayed_terms=plot_counts,
                versions=list(clusterProfiler=as.character(packageVersion('clusterProfiler')),
                  org.Hs.eg.db=as.character(packageVersion('org.Hs.eg.db'))))
record$displayed_terms <- plot_counts
record$effective_annotation_background_genes <- effective_backgrounds
record$background_interpretation <- paste(record$background_genes,'detected candidate symbols; GO/KEGG hypergeometric tests use the subset mapped to each annotation resource.')
record$provisional_annotation <- 'Fibroblast predictions are source-dominated and provisional; marker enrichment is exploratory, not an AD-specific finding.'
write_json(record,file.path(out,'enrichment_receipt.json'),auto_unbox=TRUE,pretty=TRUE)
