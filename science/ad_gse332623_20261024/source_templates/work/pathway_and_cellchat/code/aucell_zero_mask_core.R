## Explicit zero exclusion around the standard seeded AUCell ranking.
## AUCell 1.24.0's keepZeroesAsNA flag does not implement this correctly.
score_aucell_zero_mask <- function(matrix,sets,membership,aucell_seed) {
  stopifnot(identical(rownames(matrix),colnames(membership)))
  set.seed(aucell_seed)
  ranking<-AUCell_buildRankings(matrix,plotStats=FALSE,splitByBlocks=TRUE,
    BPPARAM=SerialParam(RNGseed=aucell_seed),keepZeroesAsNA=FALSE,verbose=FALSE)
  cap<-ceiling(.05*nrow(matrix))
  weights<-cap-as.matrix(getRanking(ranking));weights[weights<0]<-0
  ## Preserve all expressed-gene ranks; do not rank a smaller gene universe.
  weights<-weights*as.matrix(matrix>0)
  k<-pmin(lengths(sets),cap-1L);denominator<-k*cap-k*(k+1)/2
  scores<-as.matrix(membership%*%weights)/denominator
  dimnames(scores)<-list(names(sets),colnames(matrix))
  detected<-as.matrix(membership%*%(matrix>0))
  zero<-detected==0
  stopifnot(all(is.finite(scores)),all(scores[zero]==0),all(scores>=0),all(scores<=1))
  list(scores=scores,zero_detection_pairs=sum(zero),zero_detection_nonzero_scores=sum(scores[zero]!=0))
}

run_aucell_zero_mask_chunks <- function(cfg,branch) {
  axes<-read_axes(cfg$input);gs<-prepare_sets_alias(cfg$genesets,axes$genes,cfg$symbol_alias_path)
  dir.create(branch,recursive=TRUE,showWarnings=FALSE)
  output<-file.path(branch,'scores_three_methods_zero_mask_and_alias.h5')
  write.table(gs$registry,file.path(branch,'pathway_gene_overlap.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  write.table(data.frame(pathway_id=names(gs$sets)[gs$affected]),file.path(branch,'alias_recomputed_pathways.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  nc<-length(axes$cells);ns<-length(gs$sets);starts<-seq.int(1L,nc,by=cfg$chunk_size)
  create<-!file.exists(output);out<-H5File$new(output,if(create)'w' else 'r+')
  on.exit(out$close_all(),add=TRUE)
  if(create) {
    out[['cell_ids']]<-axes$cells;out[['pathway_ids']]<-names(gs$sets)
    h5attr(out,'source_input')<-normalizePath(cfg$input)
    h5attr(out,'zero_policy')<-'Explicit RNA>0 mask; seeded standard ranking; unchanged expressed ranks and 5% cap'
    h5attr(out,'seed')<-as.integer(cfg$seed);h5attr(out,'chunk_size')<-as.integer(cfg$chunk_size)
    sg<-out$create_group('scores')
    for(method in c('AUCell','UCell','singscore'))sg$create_dataset(method,dtype=h5types$H5T_IEEE_F32LE,
      dims=c(ns,nc),chunk_dims=c(ns,min(cfg$chunk_size,nc)),gzip_level=4)
    out[['completed_chunks']]<-integer(length(starts))
    out[['zero_detection_pairs']]<-numeric(length(starts))
    out[['zero_detection_nonzero_scores']]<-integer(length(starts));out$flush()
  } else {
    stopifnot(identical(as.character(out[['cell_ids']][]),axes$cells),identical(as.character(out[['pathway_ids']][]),names(gs$sets)))
    stopifnot(h5attr(out,'seed')==cfg$seed,h5attr(out,'chunk_size')==cfg$chunk_size)
  }
  done<-out[['completed_chunks']][];stopifnot(identical(which(done==1L),seq_len(sum(done))))
  pending<-which(done==0L)
  if(length(pending)) {
    stopifnot(cfg$workers>=1L,cfg$workers<=8L)
    cl<-parallel::makePSOCKcluster(cfg$workers,outfile='');on.exit(parallel::stopCluster(cl),add=TRUE)
    worker_pids<-unlist(parallel::clusterCall(cl,Sys.getpid))
    write_json(list(parent_R_pid=Sys.getpid(),worker_pids=worker_pids,score_RSS_budget_GiB=64,workers=cfg$workers),file.path(branch,'score_processes.json'),auto_unbox=FALSE,pretty=TRUE)
    parallel::clusterCall(cl,function(core,fix,alias_core,input,genesets,alias_path,baseline) {
      source(core,local=.GlobalEnv);source(fix,local=.GlobalEnv);source(alias_core,local=.GlobalEnv)
      assign('.worker_axes',read_axes(input),envir=.GlobalEnv)
      assign('.worker_sets',prepare_sets_alias(genesets,.worker_axes$genes,alias_path),envir=.GlobalEnv)
      assign('.worker_input',H5File$new(input,mode='r'),envir=.GlobalEnv)
      assign('.worker_baseline',H5File$new(baseline,mode='r'),envir=.GlobalEnv);NULL
    },file.path(cfg$code_dir,'pathway_core.R'),file.path(cfg$code_dir,'aucell_zero_mask_core.R'),
      file.path(cfg$code_dir,'symbol_alias_view.R'),cfg$input,cfg$genesets,cfg$symbol_alias_path,cfg$raw_scores)
    batches<-split(pending,ceiling(seq_along(pending)/cfg$workers))
    for(batch in batches) {
      while(!is.null(cfg$resource_pause_file)&&file.exists(cfg$resource_pause_file))Sys.sleep(5)
      tasks<-lapply(batch,function(k)list(k=k,start=starts[k],end=min(starts[k]+cfg$chunk_size-1L,nc),seed=as.integer(cfg$seed+k-1L)))
      result<-parallel::parLapply(cl,tasks,function(task) {
        clock<-proc.time()[['elapsed']]
        m<-read_rna_chunk(.worker_input,.worker_axes$genes,.worker_axes$cells,task$start,task$end)
        value<-score_aucell_zero_mask(m,.worker_sets$sets,.worker_sets$membership,task$seed)
        idx<-.worker_sets$affected
        corrected<-score_ucell_singscore_subset(m,.worker_sets$sets[idx],.worker_sets$membership[idx,,drop=FALSE])
        start<-task$start;end<-task$end
        for(method in c('UCell','singscore')) {
          baseline<-.worker_baseline[[paste0('scores/',method)]][,start:end]
          baseline[idx,]<-corrected[[method]];value[[method]]<-baseline
        }
        rss<-current_rss_gib();if(rss>6)stop('Corrected score worker RSS guard exceeded 6GiB')
        list(value=value,task=task,seconds=proc.time()[['elapsed']]-clock,RSS_GiB=rss)
      })
      for(item in result) {
        t<-item$task;k<-t$k;start<-t$start;end<-t$end
        out[['scores/AUCell']][,start:end]<-item$value$scores
        for(method in c('UCell','singscore'))out[[paste0('scores/',method)]][,start:end]<-item$value[[method]]
        out[['zero_detection_pairs']][k]<-item$value$zero_detection_pairs
        out[['zero_detection_nonzero_scores']][k]<-item$value$zero_detection_nonzero_scores
        out$flush();out[['completed_chunks']][k]<-1L;out$flush()
        cat(toJSON(list(branch='RNA_scores_zero_mask_and_alias_v1',chunk=t$k,chunks=length(starts),start=t$start,end=t$end,
          seconds=item$seconds,zero_detection_pairs=item$value$zero_detection_pairs,
          zero_detection_nonzero_scores=item$value$zero_detection_nonzero_scores,worker_RSS_GiB=item$RSS_GiB),auto_unbox=TRUE),'\n');flush.console()
      }
      rm(result);gc(FALSE)
    }
    parallel::clusterCall(cl,function(){.worker_input$close_all();.worker_baseline$close_all();NULL})
  }
  stopifnot(all(out[['completed_chunks']][]==1L),all(out[['zero_detection_nonzero_scores']][]==0L))
  h5attr(out,'status')<-'all_chunks_complete';out$flush()
  list(status='all_chunks_complete',branch='RNA_scores_zero_mask_and_alias_v1',cells=nc,pathways=ns,
    chunks=length(starts),zero_detection_pairs=sum(out[['zero_detection_pairs']][]),
    zero_detection_nonzero_scores=sum(out[['zero_detection_nonzero_scores']][]),
    auCell_maxRank=ceiling(.05*length(axes$genes)),seed=cfg$seed,
    zero_policy='RNA zero genes contribute zero; expressed-gene ranks, full universe and denominator unchanged',
    symbol_alias_affected_pathways=names(gs$sets)[gs$affected],
    symbol_alias_recomputed_methods=c('UCell','singscore'),symbol_alias_recomputed_columns=length(gs$affected),
    symbol_alias_unchanged_columns=ns-length(gs$affected),
    software_flag_workaround='AUCell 1.24.0 keepZeroesAsNA=TRUE indexing bug; explicit RNA>0 weight mask')
}
