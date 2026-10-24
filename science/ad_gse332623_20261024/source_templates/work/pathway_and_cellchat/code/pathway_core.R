## Exact per-cell scoring on one fixed RNA gene universe; no Seurat copy.
suppressPackageStartupMessages({library(Matrix);library(hdf5r);library(AUCell);library(UCell);library(singscore);library(BiocParallel);library(jsonlite)})
options(warn=1)
if(requireNamespace('data.table',quietly=TRUE))data.table::setDTthreads(1L)

current_rss_gib <- function() {
  text<-readLines('/proc/self/status');line<-grep('^VmRSS:',text,value=TRUE)
  as.numeric(gsub('[^0-9]','',line))/1024^2
}

read_axes <- function(path) {
  f<-H5File$new(path,mode='r');on.exit(f$close_all())
  read_index<-function(group) {
    index<-h5attr(group,'_index');as.character(group[[index]][])
  }
  genes<-read_index(f[['var']]);cells<-read_index(f[['obs']])
  stopifnot(!anyDuplicated(genes),!anyDuplicated(cells),all(nzchar(genes)),all(nzchar(cells)))
  stopifnot(h5attr(f[['X']],'encoding-type')=='csr_matrix')
  stopifnot(all(as.numeric(h5attr(f[['X']],'shape'))==c(length(cells),length(genes))))
  list(genes=genes,cells=cells)
}

read_rna_chunk <- function(handle,genes,cells,start,end) {
  ## AnnData cell-row CSR vectors are directly a genes-by-cells R CSC matrix.
  ptr<-as.numeric(handle[['X/indptr']][start:(end+1L)])
  bounds<-c(ptr[1],ptr[length(ptr)]);nnz<-diff(bounds)
  stopifnot(nnz<.Machine$integer.max,all(diff(ptr)>=0))
  if(nnz) {
    positions<-seq(from=bounds[1]+1,to=bounds[2])
    values<-as.numeric(handle[['X/data']][positions]);indices<-as.integer(handle[['X/indices']][positions])
  } else {values<-numeric();indices<-integer()}
  stopifnot(all(is.finite(values)),all(values>=0),all(indices>=0L),all(indices<length(genes)))
  new('dgCMatrix',x=values,i=indices,p=as.integer(ptr-ptr[1]),
      Dim=as.integer(c(length(genes),end-start+1L)),Dimnames=list(genes,cells[start:end]))
}

prepare_sets <- function(path,genes) {
  original<-readRDS(path);stopifnot(is.list(original),length(original)==359L)
  ids<-sub('^(hsa[0-9]{5}).*$','\\1',names(original))
  stopifnot(!anyDuplicated(ids),all(grepl('^hsa[0-9]{5}$',ids)))
  overlap<-lapply(original,function(x)intersect(unique(as.character(x)),genes))
  n<-lengths(overlap);keep<-n>=1L & n<=500L
  registry<-data.frame(pathway_id=ids,source_gene_set_name=names(original),
    pathway_name=sub('^hsa[0-9]{5}_','',names(original)),source_gene_count=lengths(original),
    RNA_gene_overlap=n,included=keep,reason=ifelse(n==0,'No RNA gene overlap',ifelse(n>500,'RNA overlap exceeds maxGSSize=500','Included')),
    stringsAsFactors=FALSE)
  selected<-overlap[keep];names(selected)<-ids[keep]
  stopifnot(length(selected)>0)
  i<-rep(seq_along(selected),lengths(selected));j<-match(unlist(selected,use.names=FALSE),genes)
  membership<-sparseMatrix(i=i,j=j,x=1,dims=c(length(selected),length(genes)),dimnames=list(names(selected),genes))
  list(sets=selected,registry=registry,membership=membership)
}

score_chunk <- function(matrix,sets,membership,aucell_seed) {
  stopifnot(identical(rownames(matrix),colnames(membership)))
  set.seed(aucell_seed)
  rank_a<-AUCell_buildRankings(matrix,plotStats=FALSE,splitByBlocks=TRUE,
    BPPARAM=SerialParam(RNGseed=aucell_seed),keepZeroesAsNA=FALSE,verbose=FALSE)
  maxrank<-ceiling(0.05*nrow(matrix))
  ## Algebraically identical to AUCell .AUC.geneSet_norm/.auc:
  ## sum(max(aucThreshold-rank,0)) / maximum possible recovery AUC.
  weights<-maxrank-as.matrix(getRanking(rank_a));weights[weights<0]<-0
  k<-pmin(lengths(sets),maxrank-1L);denominator<-k*maxrank-k*(k+1)/2
  a<-as.matrix(membership%*%weights)/denominator
  dimnames(a)<-list(names(sets),colnames(matrix));rm(rank_a,weights);gc(FALSE)
  ## Exact UCell 2.6.2 positive-signature u_stat: average descending ranks,
  ## ranks > maxRank mapped to maxRank+1, all-outside signatures score zero.
  stopifnot(!any(grepl('[+-]$',unlist(sets,use.names=FALSE))))
  max_u<-min(1500L,nrow(matrix));n_set<-lengths(sets)
  ranks_u<-matrixStats::colRanks(-as.matrix(matrix),ties.method='average',preserveShape=TRUE)
  has_top_gene<-as.matrix(membership%*%(ranks_u<=max_u))>0
  ranks_u[ranks_u>max_u]<-max_u+1L
  u<-1-(as.matrix(membership%*%ranks_u)-n_set*(n_set+1)/2)/(n_set*max_u)
  u[!has_top_gene|u<0]<-0
  dimnames(u)<-list(names(sets),colnames(matrix));rm(ranks_u,has_top_gene);gc(FALSE)
  ranks<-singscore::rankGenes(as.matrix(matrix),tiesMethod='min')
  mean_ranks<-as.matrix(membership%*%ranks)/n_set
  ## The official simpleScore(centerScore=FALSE, knownDirection=TRUE)
  ## normalization, omitting unused per-cell MAD dispersion calculations.
  low<-(n_set+1)/2
  s<-(mean_ranks-low)/(nrow(matrix)-n_set)
  dimnames(s)<-list(names(sets),colnames(matrix));rm(ranks,mean_ranks);gc(FALSE)
  stopifnot(all(is.finite(a)),all(is.finite(u)),all(is.finite(s)))
  list(AUCell=a,UCell=u,singscore=s)
}

run_pathway_chunks <- function(input,output,genesets,chunk_size=1000L,seed=2026L,checkpoints,workers=1L,code_dir=NULL,pause_file=NULL) {
  stopifnot(workers>=1L,workers<=4L)
  axes<-read_axes(input);gs<-prepare_sets(genesets,axes$genes)
  dir.create(dirname(output),recursive=TRUE,showWarnings=FALSE)
  dir.create(checkpoints,recursive=TRUE,showWarnings=FALSE)
  write.table(gs$registry,file.path(dirname(output),'pathway_gene_overlap.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  nc<-length(axes$cells);ns<-length(gs$sets);starts<-seq.int(1L,nc,by=chunk_size)
  create<-!file.exists(output);out<-H5File$new(output,if(create)'w' else 'r+')
  on.exit(out$close_all(),add=TRUE)
  if(create) {
    out[['cell_ids']]<-axes$cells;out[['pathway_ids']]<-names(gs$sets)
    h5attr(out,'source_input')<-normalizePath(input);h5attr(out,'seed')<-as.integer(seed)
    h5attr(out,'chunk_size')<-as.integer(chunk_size);h5attr(out,'RNA_gene_universe')<-as.integer(length(axes$genes))
    scoregroup<-out$create_group('scores')
    for(method in c('AUCell','UCell','singscore'))scoregroup$create_dataset(method,dtype=h5types$H5T_IEEE_F32LE,
      dims=c(ns,nc),chunk_dims=c(ns,min(chunk_size,nc)),gzip_level=4)
    out[['completed_chunks']]<-integer(length(starts));out$flush()
  } else {
    stopifnot(identical(as.character(out[['cell_ids']][]),axes$cells),identical(as.character(out[['pathway_ids']][]),names(gs$sets)))
    stopifnot(h5attr(out,'chunk_size')==chunk_size,h5attr(out,'seed')==seed)
  }
  done<-out[['completed_chunks']][];stopifnot(identical(which(done==1L),seq_len(sum(done))))
  pending<-which(done==0L)
  if(length(pending)) {
    cl<-parallel::makePSOCKcluster(workers,outfile='')
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(core,input,genesets) {
      source(core,local=.GlobalEnv)
      assign('.worker_axes',read_axes(input),envir=.GlobalEnv)
      assign('.worker_sets',prepare_sets(genesets,.worker_axes$genes),envir=.GlobalEnv)
      assign('.worker_input',H5File$new(input,mode='r'),envir=.GlobalEnv)
      NULL
    },file.path(code_dir,'pathway_core.R'),input,genesets)
    batches<-split(pending,ceiling(seq_along(pending)/workers))
    for(batch in batches) {
      while(!is.null(pause_file) && file.exists(pause_file))Sys.sleep(5)
      tasks<-lapply(batch,function(k)list(k=k,start=starts[k],end=min(starts[k]+chunk_size-1L,nc),seed=as.integer(seed+k-1L)))
      result<-parallel::parLapply(cl,tasks,function(task) {
        clock<-proc.time()[['elapsed']]
        m<-read_rna_chunk(.worker_input,.worker_axes$genes,.worker_axes$cells,task$start,task$end)
        scores<-score_chunk(m,.worker_sets$sets,.worker_sets$membership,task$seed)
        rss<-current_rss_gib();if(rss>20)stop('Worker RSS guard exceeded 20GiB')
        list(scores=scores,info=list(chunk=task$k,start=task$start,end=task$end,seconds=proc.time()[['elapsed']]-clock,worker_RSS_GiB=rss),seed=task$seed)
      })
      for(item in result) {
        k<-item$info$chunk;start<-item$info$start;end<-item$info$end
        for(method in c('AUCell','UCell','singscore'))out[[paste0('scores/',method)]][,start:end]<-item$scores[[method]]
        out$flush()
        saveRDS(list(chunk=k,start=start,end=end,aucell_seed=item$seed),file.path(checkpoints,sprintf('chunk_%05d.rds',k)))
        out[['completed_chunks']][k]<-1L;out$flush()
        info<-c(item$info,list(chunks=length(starts),workers=workers,parent_RSS_GiB=current_rss_gib()))
        cat(toJSON(info,auto_unbox=TRUE),'\n');flush.console()
      }
      rm(result);gc(FALSE)
    }
    parallel::clusterCall(cl,function(){.worker_input$close_all();NULL})
  }
  h5attr(out,'status')<-'all_chunks_complete';out$flush()
  list(status='all_chunks_complete',cells=nc,pathways=ns,genes=length(axes$genes),chunks=length(starts),
    workers=workers,methods=c('AUCell','UCell','singscore'),seed=seed,AUCell_chunk_seed='seed + chunk_index - 1; explicit SerialParam RNGseed',aucellMaxRank=ceiling(0.05*length(axes$genes)),
    ucellMaxRank=min(1500,length(axes$genes)),singscore_centered=FALSE,
    excluded_gene_sets=gs$registry[!gs$registry$included,],RSS_GiB=current_rss_gib())
}
