args<-commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1L)
cfg<-jsonlite::fromJSON(args[1],simplifyVector=TRUE)
source(file.path(cfg$code_dir,'pathway_core.R'))
result<-run_pathway_chunks(cfg$input,cfg$raw_scores,cfg$genesets,as.integer(cfg$chunk_size),2026L,cfg$checkpoints,as.integer(cfg$workers),cfg$code_dir,cfg$resource_pause_file)
result$package_versions<-as.list(vapply(c('AUCell','UCell','singscore','Matrix','hdf5r','BiocParallel'),function(x)as.character(packageVersion(x)),character(1)))
result$R_version<-R.version.string
jsonlite::write_json(result,file.path(cfg$output,'r_scoring_completion.json'),auto_unbox=TRUE,pretty=TRUE)
