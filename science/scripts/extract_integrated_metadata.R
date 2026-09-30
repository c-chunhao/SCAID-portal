suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(jsonlite)})
b <- Sys.getenv('SCAID_SCIENCE_TABLES')
for(nm in c('B_Plasma','T_NK')) {
 p <- if(nm=='T_NK') Sys.getenv('SCAID_TNK_INTEGRATED_RDS') else Sys.getenv('SCAID_BPLASMA_INTEGRATED_RDS')
 message(format(Sys.time()),' Reading ',nm,' original integration object')
 x<-readRDS(p);md<-x@meta.data
 cols<-intersect(c('OldCells','Sample','Disease','DiseaseFull','Dataset','SampleType','SampleTypeFull','MajorCellType','celltype','celltype1','celltype2'),names(md));md<-md[,cols,drop=FALSE];md$integration_cell_id<-rownames(md)
 write.table(md,gzfile(file.path(b,paste0(nm,'_original_integration_metadata.tsv.gz'))),sep='\t',quote=FALSE,row.names=FALSE,na='')
 write_json(list(source=p,n_cells=nrow(md),columns=names(md),note='Metadata-only extraction from original published integration source; counts and annotation unchanged.'),file.path(b,paste0(nm,'_metadata_evidence.json')),auto_unbox=TRUE,pretty=TRUE)
 message(format(Sys.time()),' Saved ',nm,' metadata ',nrow(md));rm(x,md);gc()
}
