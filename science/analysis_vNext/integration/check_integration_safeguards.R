suppressPackageStartupMessages({library(SeuratObject);library(Matrix);library(jsonlite)})
script_file<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
base<-dirname(normalizePath(script_file,mustWork=TRUE))
source(file.path(base,'integration_safeguards.R'))
freeze_analysis_rng()
c <- Matrix(matrix(c(1,0,3,2,1,0,0,4,2,1,0,1),nrow=3),sparse=TRUE)
rownames(c)<-c('CD3D','MS4A1','LYZ');colnames(c)<-paste0('cell',1:4)
x <- CreateSeuratObject(c,min.cells=0,min.features=0)
retrieved <- get_RNA_counts_checked(x)
stopifnot(isTRUE(all.equal(retrieved,c,check.attributes=FALSE)))
emb<-matrix(seq_len(12),nrow=4,dimnames=list(rev(colnames(c)),paste0('PC_',1:3)))
x<-attach_reduction_checked(x,emb,'aligned')
stopifnot(identical(rownames(Embeddings(x,'aligned')),Cells(x)),
          isTRUE(all.equal(unname(Embeddings(x,'aligned')),unname(emb[match(Cells(x),rownames(emb)),,drop=FALSE]),check.attributes=FALSE)))
bad <- emb;rownames(bad)[1] <- 'absent'
rejected_missing <- inherits(try(attach_reduction_checked(x,bad,'bad'),silent=TRUE),'try-error')
md1<-data.frame(celltype1=c('Tn','Te'),celltype2=c('Tn_1','Te_1'),row.names=c('cell2','cell1'))
md2<-data.frame(celltype1=c('Bn','Mono'),celltype2=c('Bn_1','Mono_1'),row.names=c('cell4','cell3'))
a<-assemble_annotations_checked(x,list(md1,md2))
md<-rbind(md1,md2)
stopifnot(identical(as.character(a$celltype1),as.character(md[Cells(a),'celltype1'])))
rejected_duplicate <- inherits(try(assemble_annotations_checked(x,list(md1,md1)),silent=TRUE),'try-error')
norm<-full_RNA_log_normalization(c)
stopifnot(max(abs(Matrix::colSums(expm1(norm))-10000))<1e-6)
# The denominator remains the complete library after a feature subset.
subset_one_gene <- norm[1,,drop=FALSE]
stopifnot(all(Matrix::colSums(expm1(subset_one_gene))<=10000),any(Matrix::colSums(expm1(subset_one_gene))<10000))
stopifnot(rejected_missing,rejected_duplicate)
# Store scaled features in an order different from the raw features, so the
# exporter cannot silently substitute a raw-gene membership-filtered order.
x<-SetAssayData(x,assay='RNA',layer='scale.data',new.data=as.matrix(norm[c(3,1),,drop=FALSE]))
tmp_source<-tempfile(fileext='.rds');saveRDS(x,tmp_source)
tmp_axes<-tempfile(fileext='.json')
export_RNA_layer_axes_checked(x,tmp_axes,tmp_source)
axes<-read_json(tmp_axes,simplifyVector=TRUE)
actual_scale_axis<-rownames(GetAssayData(x,assay='RNA',layer='scale.data'))
stopifnot(identical(axes$scaled_HVG_feature_names,actual_scale_axis),
          identical(axes$raw_feature_names,rownames(get_RNA_counts_checked(x))),
          identical(sort(axes$scaled_HVG_feature_names),sort(c('LYZ','CD3D'))))
unlink(c(tmp_source,tmp_axes))
write_json(list(reduction_permutation_reordered_by_ID=TRUE,missing_reduction_cells_rejected=rejected_missing,
                duplicate_annotation_cells_rejected=rejected_duplicate,annotation_metadata_reordered_by_ID=TRUE,
                full_sparse_library_normalizes_to_10000=TRUE,HVG_subset_does_not_renormalize=TRUE,
                SeuratObject_layer_API_counts_retrieval=TRUE,
                exported_scale_feature_axis_keeps_exact_source_row_order=TRUE),
           file.path(base,'integration_safeguards_validation.json'),auto_unbox=TRUE,pretty=TRUE)
cat('Scientific cell-axis and normalization guards passed\n')
