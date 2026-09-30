suppressPackageStartupMessages(library(jsonlite))
p <- Sys.getenv('SCAID_RELEASE_METADATA_RDS')
b <- Sys.getenv('SCAID_SCIENCE_TABLES')
x<-readRDS(p);print(names(x))
z<-lapply(seq_along(x),function(i){m<-x[[i]]; m<-m[,intersect(c('OldCells','Sample','Disease','DiseaseFull','Dataset','SampleType','SampleTypeFull','MajorCellType'),names(m)),drop=FALSE];m$Data_object<-names(x)[i];m})
m<-do.call(rbind,z);rownames(m)<-NULL
write.table(m,gzfile(file.path(b,'all_release_original_metadata.tsv.gz')),sep='\t',quote=FALSE,row.names=FALSE,na='')
write_json(list(source=p,n_cells=nrow(m),n_objects=length(x),columns=names(m),object_cell_counts=sapply(x,nrow),note='Original complete released object metadata, legacy filenames/condition labels retained until source correction in Python.'),file.path(b,'all_release_metadata_evidence.json'),auto_unbox=TRUE,pretty=TRUE)
