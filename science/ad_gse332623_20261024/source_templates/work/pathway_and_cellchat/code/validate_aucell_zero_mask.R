args<-commandArgs(TRUE);stopifnot(length(args)==4L)
source(file.path(args[1],'pathway_core.R'));source(file.path(args[1],'aucell_zero_mask_core.R'))
axes<-read_axes(args[2]);gs<-prepare_sets('@PROJECT@/20240526_kegg_hsa_list.rds',axes$genes)
fh<-H5File$new(args[2],'r');m<-read_rna_chunk(fh,axes$genes,axes$cells,1L,1000L);fh$close_all()
fast<-score_aucell_zero_mask(m,gs$sets,gs$membership,2026L)
reference<-H5File$new(args[3],'r');stopifnot(identical(as.character(reference[['cell_id']][]),colnames(m)),identical(as.character(reference[['pathway_id']][]),names(gs$sets)))
official<-reference[['explicit_zero_mask']][,];reference$close_all()
error<-max(abs(fast$scores-official));stopifnot(error<1e-12,fast$zero_detection_nonzero_scores==0L)
answer<-list(status='validated',cells=1000L,pathways=length(gs$sets),seed=2026L,aucMaxRank=1683L,
  explicit_RNA_zero_to_NA_official_calcAUC_max_absolute_error=error,
  zero_detection_pairs=fast$zero_detection_pairs,zero_detection_nonzero_scores=fast$zero_detection_nonzero_scores,
  reference=args[3],reference_semantics='Independent official AUCell ranking with RNA-zero positions set NA before official calcAUC',
  zero_policy='Seeded standard full-gene ranking weights multiplied by RNA>0; expressed-gene rank and denominator unchanged')
write_json(answer,args[4],pretty=TRUE,auto_unbox=TRUE,digits=17)
cat('INDEPENDENT_ZERO_MASK_VALIDATED\n')
