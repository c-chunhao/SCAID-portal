suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(SeuratObject))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)<3L) stop('Usage: Rscript export_final_reductions.R ATLAS OUTPUT_DIR SOURCE_RDS')
atlas <- args[1]; dest <- args[2]; src <- args[3]
stopifnot(atlas %in% c('T_NK','B_Plasma'), file.exists(src))
dir.create(dest, recursive=TRUE, showWarnings=FALSE)
cat('Loading',atlas,src,format(Sys.time(),tz='Asia/Shanghai'), '\n');flush.console()
x <- readRDS(src)
meta <- x@meta.data
expected <- if(atlas=='T_NK') 937825L else 253635L
stopifnot(nrow(meta)==expected,!anyDuplicated(rownames(meta)))
cat('Loaded metadata',nrow(meta),'reductions',names(x@reductions),'\n');flush.console()
reductions <- list()
for(nm in names(x@reductions)) {
 z <- x@reductions[[nm]]@cell.embeddings
 same <- identical(rownames(z),rownames(meta))
 item <- list(name=nm,dimensions=ncol(z),cells=nrow(z),same_order_as_metadata=same,finite=all(is.finite(z)),class=class(x@reductions[[nm]]))
 if(ncol(z)>2 && nrow(z)==nrow(meta)) {
  stopifnot(same,all(is.finite(z)))
  out <- data.frame(integration_cell_id=rownames(z),z,check.names=FALSE)
  path <- file.path(dest,paste0(atlas,'_final_',nm,'.tsv.gz'))
  write.table(out,gzfile(path),sep='\t',row.names=FALSE,quote=FALSE)
  item$export_path <- path
 }
 reductions[[nm]] <- item
}
meta_out <- meta[,intersect(c('OldCells','Sample','Disease','DiseaseFull','Dataset','SampleType','SampleTypeFull','MajorCellType','celltype','celltype1','celltype2'),names(meta)),drop=FALSE]
meta_out$integration_cell_id <- rownames(meta)
write.table(meta_out,gzfile(file.path(dest,paste0(atlas,'_final_metadata.tsv.gz'))),sep='\t',quote=FALSE,row.names=FALSE)
commands <- lapply(x@commands,function(c) list(name=c@name,time=as.character(c@time.stamp),assay=c@assay.used,call=c@call.string,params=c@params))
write_json(list(atlas=atlas,source=src,source_bytes=unname(file.info(src)$size),n_cells=nrow(meta),metadata_columns=names(meta),reductions=reductions,commands=commands,R=as.character(getRversion()),SeuratObject=as.character(packageVersion('SeuratObject')),source_unchanged=TRUE,no_reintegration=TRUE),file.path(dest,paste0(atlas,'_final_object_inspection.json')),pretty=TRUE,auto_unbox=TRUE,null='null')
cat('Finished',atlas,format(Sys.time(),tz='Asia/Shanghai'),'\n')
