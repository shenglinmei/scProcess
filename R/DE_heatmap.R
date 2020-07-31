
DEheatmap=function(p2myeloid,conosCluster,appname,removeGene=NULL,num=20){
  de1 <- p2myeloid$getDifferentialGenes(groups=conosCluster)
  ggg=NULL
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
      index=grepl('^MT-',rownames(z))
      z=z[!index,]
      index=grepl('^HSP',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    #markers=markers[!is.na(markers)]
    
    #allg=cbind(allg,markers)
    ggg=c(ggg,markers)
    
  }
  
  markers=unique(ggg)
  
  
  
  
  
  x <- as.matrix(p2myeloid$counts[names(conosCluster),markers])
  ## trim outliers
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.01,0.98))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  ## sample 2000 cells for plotting
  #x <- x[,sample(ncol(x),2000)]
  o <- order(conosCluster[colnames(x)])
  
  
  x=x[,o]
  ss=apply(x,2,function(x) sum(x))
  
  x=x[,abs(ss)>0.1]
  
  
  annot <- data.frame(CellType=conosCluster[colnames(x)],row.names = colnames(x))
  
  
  annot2=annot
  
  #annot2 = data.frame(ID = as.factor(as.character(annot[,1])))
  rownames(annot2)=colnames(x)
  
  
  cellA=annot2[,1]
  names(cellA)=rownames(annot2)
  
  o <- order(cellA)
  
  x=x[,o]
  
  pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
  ## draw heatmap
  
  fout=paste(appname,'.DiffG.png')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  heat=pheatmap(x,cluster_cols=FALSE,annotation_col = annot2,show_colnames = F,annotation_legend = TRUE, #,gaps_col =  1,
                cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =5,width=8,height=8.8,   #3.5*0.02*length(markers),
                breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
  
}




DEcaculate2.term=function(p2,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  allg=NULL
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  folderN=paste(appname,'.DE',sep='')
  pwd=getwd()
  pwd2=paste(pwd,'/',folderN,'/',sep='')
  
  system(paste('mkdir ',folderN))
  
  setwd(pwd2)
  
  
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      
      figGO=GOanalysis.term(markers,iid)
      figGO=figGO[['BP']]
      
      for (jj in seq(nrow(figGO))){
        
        gs=figGO[jj,'genes']
        gs=strsplit(gs,',')[[1]]
        
        if (length(gs) >= 4 ){
          
          fout=paste(appname,'_',n,'_',figGO[jj,'Term'],'_BP.heatmap.png',sep='')
          rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
          pheatmap(x[gs,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
                   cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=5*0.02*length(gs),
                   breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
          
        }}
      
      
    }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  setwd(pwd)
}



DEcaculate2=function(p2,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  allg=NULL  
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.new.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      GOanalysis(markers,iid) }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  
}  






DEheatmap=function(con,groups,anoSample,prefix,num=20,z.threshold=3,pagodaDE=NULL){
  # con: conos object
  # anoSample: sample annotation
  # prefix: prefix of output figure
  # num: number of higly expressed gene
  # pagodaDE: for pagoda object
  
  # get DE genes
  if (!is.null(pagodaDE)){
    cds=mergeMatrix(lapply(con$samples,function(x) t(x$misc$rawCounts)))
    tmp2=p2 <- Pagoda2$new(cds)
    tmp2$adjustVariance(plot = T, gam.k = 10)
    de1 <- tmp2$getDifferentialGenes(groups=groups,z.threshold=z.threshold)
  }else{
    de1 <- con$getDifferentialGenes(groups=groups,z.threshold=z.threshold)
    
  }
  
  # top num higly expressed genes
  gs=unlist(lapply(de1,function(x) { x=x[order(x$Z,decreasing=T),]
  rownames(x)[1:min(c(nrow(x),num))]
  }))
  
  gs=gs[!is.na(gs)]
  
  
  # get expression matrix
  exp=mergeMatrix(lapply(con$samples,function(x) t(x$counts)))
  
  
  x <- t(as.matrix(exp[gs,names(groups)]))
  ## trim outliers
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.01,0.98))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  
  # transform to Z-score
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  
  o <- order(groups[colnames(x)])
  
  x=x[,o]
  annot <- data.frame(CellType=groups[colnames(x)],'Sample'=anoSample[colnames(x)],row.names = colnames(x))
  
  pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
  ## draw heatmap
  
  fout=paste(prefix,'.DiffG.png',sep='')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  heat=pheatmap(x,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,annotation_legend = TRUE, #,gaps_col =  1,
                cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =5,width=8,height=3.5*0.02*length(gs),
                breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
  
}





DEcaculate2.conos=function(con,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  
  if (!requireNamespace("pagoda2", quietly = TRUE)) {
    stop("You have to install pagoda2")
  }
  allg=NULL
  
  cds=mergeMatrix(lapply(con$samples,function(x) t(x$misc$rawCounts)))
  cds=cds[,names(conosCluster)]
  p2 <- Pagoda2$new(cds)
  p2$adjustVariance(plot = T, gam.k = 10)
  
  de1 <- p2$getDifferentialGenes(groups=conosCluster,z.threshold = cutoff)
  
  
  f1=paste(appname,'_diffGene.rds',sep='')
  saveRDS(de1,f1)
  for(n in names(de1)){
    x=de1[[n]]
    z <- x[order(-x$Z),]
    if(!is.null(removeGene)){
      index=grepl('^RP[LKS]',rownames(z))
      z=z[!index,]
    }
    markers=rownames(z)[1:num]
    markers=markers[!is.na(markers)]
    allg=cbind(allg,markers)
    x <- as.matrix(p2$counts[names(conosCluster),markers])
    ## trim outliers
    x <- apply(x, 2, function(xp) {
      qs <- quantile(xp,c(0.01,0.98))
      xp[xp<qs[1]] <- qs[1]
      xp[xp>qs[2]] <- qs[2]
      xp
    })
    x <- x[,(apply(x,2,sd) != 0)]
    x <- t(scale(x))
    ## sample 2000 cells for plotting
    #x <- x[,sample(ncol(x),2000)]
    o <- order(conosCluster[colnames(x)])
    annot <- data.frame(cluster=conosCluster[colnames(x)],row.names = colnames(x))
    
    annot$dtype='other'
    annot[as.character(annot[,1])==n,'dtype']=n
    annot$dtype=as.factor(annot$dtype)
    
    
    pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
    ## draw heatmap
    
    fout=paste(appname,'_',n,'_marker.heatmap.new.png',sep='')
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    pheatmap(x[,o],labels_col=FALSE,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=5,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      GOanalysis(markers,iid) }
  }
  colnames(allg)=names(de1)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  
}


#DEheatmap(con,groups,anoSample,'conosDE',num=20,z.threshold=0)



