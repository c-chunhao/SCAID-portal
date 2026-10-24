render_cellchat<-function(cellchat,out,cfg) {
  set.seed(cfg$seed)
  figure_dir<-file.path(out,'figures');dir.create(figure_dir,showWarnings=FALSE)
  groupSize<-as.numeric(table(cellchat@idents));names(groupSize)<-names(table(cellchat@idents))
  manifest<-list();errors<-list();idx<-0L;erridx<-0L
  role_data<-list();role_index<-0L
  role_scatter<-function(signaling=NULL,title) {
    p<-netAnalysis_signalingRole_scatter(cellchat,signaling=signaling,title=title,label.size=3,font.size=11)
    prob<-cellchat@net$prob;pval<-cellchat@net$pval
    if(!is.null(signaling)) {
      ids<-dimnames(prob)[[3]]
      idx<-which(cellchat@DB$interaction[ids,'pathway_name']%in%signaling)
      prob<-prob[,,idx,drop=FALSE];pval<-pval[,,idx,drop=FALSE]
    }
    count<-apply(prob>0 & pval<0.05,c(1,2),sum)
    links<-rowSums(count)+colSums(count)-diag(count)
    p$data$Count<-as.numeric(links[as.character(p$data$labels)])
    stopifnot(!anyNA(p$data$Count),all(p$data$Count>=0))
    if(any(prob>0 & pval<0.05))stopifnot(any(p$data$Count>0))
    role_index<<-role_index+1L
    role_data[[role_index]]<<-data.frame(signaling=if(is.null(signaling))'All' else signaling,cell_type=as.character(p$data$labels),incoming_strength=p$data$y,outgoing_strength=p$data$x,LR_count=p$data$Count)
    p+theme(plot.margin=margin(16,24,16,16))
  }
  safe<-function(s)gsub('[^A-Za-z0-9._-]','_',s)
  plot_file<-function(family,id,title,fun,width=8,height=8) {
    target_dir<-file.path(figure_dir,family);dir.create(target_dir,recursive=TRUE,showWarnings=FALSE)
    path<-file.path(target_dir,paste0(safe(id),'.png'))
    tryCatch({
      png(path,width=width,height=height,units='in',res=180,type='cairo',bg='white');on.exit(if(dev.cur()>1)dev.off(),add=TRUE)
      par(mar=c(2,2,4,2),xpd=TRUE);fun();dev.off()
      idx<<-idx+1L;manifest[[idx]]<<-data.frame(object_id=cfg$object_id,family=family,figure_id=id,title=title,path=path,width_px=as.integer(width*180),height_px=as.integer(height*180),stringsAsFactors=FALSE)
    },error=function(e){
      if(dev.cur()>1)dev.off();if(file.exists(path))unlink(path)
      erridx<<-erridx+1L;errors[[erridx]]<<-data.frame(family=family,figure_id=id,error=conditionMessage(e))
    })
  }
  plot_file('Overview','Count','Number of interactions',function()netVisual_circle(cellchat@net$count,vertex.weight=groupSize,weight.scale=TRUE,label.edge=FALSE,title.name='Number of interactions',vertex.label.cex=.85,margin=.35))
  plot_file('Overview','Weight','Interaction strength',function()netVisual_circle(cellchat@net$weight,vertex.weight=groupSize,weight.scale=TRUE,label.edge=FALSE,title.name='Interaction strength',vertex.label.cex=.85,margin=.35))
  plot_file('Overview','cellRole','Outgoing and incoming interaction strength',function()print(role_scatter(title='All signaling pathways')),9,7)
  for(cel in levels(cellchat@idents))for(measure in c('Count','Weight')) {
    mat<-cellchat@net[[tolower(measure)]];single<-mat*0;single[cel,]<-mat[cel,]
    title<-paste(cel,if(measure=='Count')'interaction counts' else 'interaction strength')
    plot_file('Single',paste0(measure,'-',cel),title,function()netVisual_circle(single,vertex.weight=groupSize,weight.scale=TRUE,edge.weight.max=max(mat),vertex.label.cex=.85,title.name=title,margin=.35))
  }
  df<-subsetCommunication(cellchat);pathways<-unique(as.character(df$pathway_name))
  if(cfg$plot_limit>0)pathways<-head(pathways,cfg$plot_limit)
  for(pw in pathways) {
    plot_file('SignalingCellRole',pw,pw,function()print(role_scatter(signaling=pw,title=pw)),9,7)
    plot_file('SignalingCircle',pw,paste(pw,'signaling pathway network'),function()netVisual_aggregate(cellchat,signaling=pw,layout='circle',signaling.name=pw,vertex.weight=groupSize,vertex.label.cex=.85,margin=.35))
    plot_file('SignalingHt',pw,paste(pw,'interaction strength'),function()ComplexHeatmap::draw(netVisual_heatmap(cellchat,signaling=pw,color.heatmap='Reds',title.name=paste(pw,'interaction strength'),font.size=11,font.size.title=13)),9,8)
    plot_file('SignalingChord',pw,paste(pw,'signaling pathway network'),function()netVisual_aggregate(cellchat,signaling=pw,layout='chord',signaling.name=pw,small.gap=1,big.gap=8,lab.cex=.85),10,10)
    n_pairs<-length(unique(df$interaction_name[df$pathway_name==pw]))
    plot_file('LRcontribution',pw,paste(pw,'ligand-receptor contributions'),function()print(netAnalysis_contribution(cellchat,signaling=pw)+ggtitle(paste(pw,'ligand-receptor contributions'))+theme(plot.margin=margin(16,20,16,16))),8,min(14,max(4,2+.22*n_pairs)))
    plot_file('SignalingRole',pw,paste(pw,'signaling roles'),function()netAnalysis_signalingRole_network(cellchat,signaling=pw,width=8,height=2.8,font.size=10,font.size.title=12),10,4)
    pairs<-extractEnrichedLR(cellchat,signaling=pw,geneLR.return=FALSE)
    for(lr in as.character(pairs[[1]])) {
      plot_file('LRcircle',paste(pw,lr,sep='__'),paste(pw,lr,sep=': '),function()netVisual_individual(cellchat,signaling=pw,layout='circle',pairLR.use=lr,signaling.name=paste(pw,lr,sep=': '),vertex.weight=groupSize,vertex.label.cex=.85,margin=.35))
    }
  }
  manifest<-if(length(manifest))do.call(rbind,manifest) else data.frame()
  write.table(manifest,file.path(out,'cellchat_figure_manifest.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  err<-if(length(errors))do.call(rbind,errors) else data.frame(family=character(),figure_id=character(),error=character())
  write.table(err,file.path(out,'cellchat_figure_errors.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  if(length(role_data))write.table(do.call(rbind,role_data),file.path(out,'cellchat_role_scatter_data.tsv'),sep='\t',quote=FALSE,row.names=FALSE)
  jsonlite::write_json(list(status=if(nrow(err))'figures_need_review' else 'figures_generated',figures=nrow(manifest),errors=nrow(err),families=as.list(table(manifest$family)),pathways_rendered=pathways,plot_limit=cfg$plot_limit),file.path(out,'cellchat_rendering.json'),pretty=TRUE,auto_unbox=TRUE)
  if(nrow(err))stop('CellChat computation finished; some plot families require repair. See cellchat_figure_errors.tsv')
  cat('CELLCHAT_FIGURES_GENERATED\n')
}
