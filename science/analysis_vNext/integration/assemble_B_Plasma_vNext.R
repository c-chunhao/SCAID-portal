# Explicit replacement for the historical Mult_04 variable-name mismatch.
# Run only on a machine with the required memory; preserve all original RDS.
suppressPackageStartupMessages({library(Seurat);library(jsonlite)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: Rscript assemble_B_Plasma_vNext.R INPUT_ROOT DISTINCT_OUTPUT_RDS')
root<-normalizePath(args[1],mustWork=TRUE)
destination<-normalizePath(args[2],mustWork=FALSE)
if(file.exists(destination))stop('Refusing to overwrite an existing RDS output')
script_file<-sub('^--file=','',commandArgs()[grepl('^--file=',commandArgs())][1])
source(file.path(dirname(normalizePath(script_file,mustWork=TRUE)),'integration_safeguards.R'))
settings<-freeze_analysis_rng()
branch_files<-file.path(root,c('02.scAnn/Mult_03/B/scRNA.B.rds','02.scAnn/Mult_03/Plasma/scRNA.Plasma.rds','02.scAnn/Mult_03/Bgc/scRNA.Bgc.rds'))
base_file<-file.path(root,'02.scAnn/Mult_01/scRNA.Mult1.rds')
# Read one branch at a time, retain metadata only, release its expression.
metadata<-lapply(branch_files,function(p){
  z<-readRDS(p);md<-z@meta.data[,c('celltype1','celltype2'),drop=FALSE]
  rm(z);gc();md
})
base_object<-readRDS(base_file)
x<-assemble_annotations_checked(base_object,metadata)
rm(base_object,metadata);gc()
# Explicit reduction/dimensions/seed, retaining counts and legacy labels.
x<-display_umap_checked(x,reduction='harmony',dims=1:20,name='UMAP_vNext',seed=20261005L)
DefaultAssay(x)<-'RNA'
saveRDS(x,destination)
write_reproducibility_record(paste0(destination,'.provenance.json'),c(base_file,branch_files),
                            c(settings,list(reduction='harmony',dims=1:20,annotation_columns=c('celltype1','celltype2'))))
