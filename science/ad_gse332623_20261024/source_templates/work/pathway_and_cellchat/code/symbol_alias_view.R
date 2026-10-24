## A verified one-to-one internal matching view; the source RNA file is immutable.
canonical_rna_view <- function(genes,path) {
  tab<-read.delim(path,stringsAsFactors=FALSE,check.names=FALSE)
  if(!'alias_confirmed'%in%names(tab))tab$alias_confirmed<-tab$source_alias_confirmed
  if(!'ensembl_exact_match'%in%names(tab))tab$ensembl_exact_match<-vapply(seq_len(nrow(tab)),function(i)tab$source_ensembl_id[i]%in%strsplit(tab$canonical_reference_Ensembl_ids[i],';',fixed=TRUE)[[1]],logical(1))
  stopifnot(nrow(tab)>0L,!anyDuplicated(tab$source_symbol),!anyDuplicated(tab$canonical_symbol),
    all(as.logical(tab$alias_confirmed)),all(as.logical(tab$ensembl_exact_match)),
    !any(as.logical(tab$canonical_present_in_RNA)),all(tab$source_symbol%in%genes),!any(tab$canonical_symbol%in%genes))
  canonical<-genes;canonical[match(tab$source_symbol,genes)]<-tab$canonical_symbol
  stopifnot(!anyDuplicated(canonical))
  list(genes=canonical,mapping=tab)
}

prepare_sets_alias <- function(path,genes,alias_path) {
  old<-prepare_sets(path,genes);view<-canonical_rna_view(genes,alias_path)
  original<-readRDS(path)
  reverse<-setNames(view$mapping$source_symbol,view$mapping$canonical_symbol)
  ## Prioritize direct source hits, then resolve absent canonical members to
  ## verified source rows; unique() counts each physical RNA gene only once.
  overlap<-lapply(original,function(x) {
    members<-unique(as.character(x));missing<-!members%in%genes & members%in%names(reverse)
    members[missing]<-unname(reverse[members[missing]])
    intersect(unique(members),genes)
  })
  count<-lengths(overlap);keep<-count>=1L & count<=500L
  stopifnot(identical(unname(keep),old$registry$included))
  new<-old;new$sets<-overlap[keep];names(new$sets)<-old$registry$pathway_id[keep]
  new$registry$RNA_gene_overlap<-count
  i<-rep(seq_along(new$sets),lengths(new$sets));j<-match(unlist(new$sets,use.names=FALSE),genes)
  new$membership<-sparseMatrix(i=i,j=j,x=1,dims=c(length(new$sets),length(genes)),dimnames=list(names(new$sets),genes))
  stopifnot(identical(names(new$sets),names(old$sets)))
  affected<-which(vapply(seq_along(new$sets),function(i)!identical(sort(new$sets[[i]]),sort(old$sets[[i]])),logical(1)))
  stopifnot(length(affected)>0L)
  new$affected<-affected;new$alias_mapping<-view$mapping
  colnames(new$membership)<-genes
  new$registry$baseline_RNA_gene_overlap<-old$registry$RNA_gene_overlap
  new$registry$alias_restored_member_count<-new$registry$RNA_gene_overlap-old$registry$RNA_gene_overlap
  new$registry$low_RNA_gene_coverage<-new$registry$RNA_gene_overlap<10L
  new
}

score_ucell_singscore_subset <- function(matrix,sets,membership) {
  max_u<-min(1500L,nrow(matrix));sizes<-lengths(sets)
  ranks<-matrixStats::colRanks(-as.matrix(matrix),ties.method='average',preserveShape=TRUE)
  top<-as.matrix(membership%*%(ranks<=max_u))>0;ranks[ranks>max_u]<-max_u+1L
  u<-1-(as.matrix(membership%*%ranks)-sizes*(sizes+1)/2)/(sizes*max_u)
  u[!top|u<0]<-0;rm(ranks,top);gc(FALSE)
  ranks<-singscore::rankGenes(as.matrix(matrix),tiesMethod='min')
  s<-(as.matrix(membership%*%ranks)/sizes-(sizes+1)/2)/(nrow(matrix)-sizes)
  dimnames(u)<-dimnames(s)<-list(names(sets),colnames(matrix))
  stopifnot(all(is.finite(u)),all(is.finite(s)))
  list(UCell=u,singscore=s)
}
