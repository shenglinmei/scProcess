plotDEheatmap2 <- function(con,groups,de=NULL,min.auc=NULL,min.specificity=NULL,min.precision=NULL,n.genes.per.cluster=10,additional.genes=NULL,labeled.gene.subset=NULL, expression.quantile=0.99,pal=colorRampPalette(c('dodgerblue1','grey95','indianred1'))(1024),ordering='-AUC',column.metadata=NULL,show.gene.clusters=TRUE, remove.duplicates=TRUE, column.metadata.colors=NULL, show.cluster.legend=TRUE, show_heatmap_legend=FALSE, border=TRUE, return.details=FALSE, row.label.font.size=10, order.clusters=FALSE, split=FALSE, split.gap=0, ...) {
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    stop("pheatmap package needs to be installed to use plotDEheatmap")
  }
  
# 
  # groups=as.factor(typefc)
  # de=sannot.de
  # min.auc=0.75
  # min.specificity=NULL
  # min.precision=NULL
  # n.genes.per.cluster=20
  # additional.genes=NULL
  # labeled.gene.subset=NULL
  # expression.quantile=0.99
  # pal=colorRampPalette(c('dodgerblue1','grey95','indianred1'))(1024)
  # ordering='-AUC'
  # column.metadata=NULL
  # show.gene.clusters=TRUE
  # remove.duplicates=TRUE
  # column.metadata.colors=NULL
  # show.cluster.legend=TRUE
  # show_heatmap_legend=FALSE
  # border=TRUE
  # return.details=FALSE
  # row.label.font.size=3
  # order.clusters=FALSE
  # split=FALSE
  # split.gap=0
  # 
  # show.gene.clusters=T
  # column.metadata=list(samples=con$getDatasetPerCell())
  # #column.metadata.colors = list(clusters=sannot.pal,samples=sample.pal)
  
  
  if(is.null(de)) { # run DE
    de <- con$getDifferentialGenes(groups=groups,append.auc=TRUE,z.threshold=0,upregulated.only=TRUE)
  }
  
  # drop empty results
  de <- de[unlist(lapply(de,nrow))>0]
  
  # drop results that are not in the factor levels
  de <- de[names(de) %in% levels(groups)]
  
  # order de list to match groups order
  de <- de[order(match(names(de),levels(groups)))]
  
  print(names(de))
  # apply filters
  if(!is.null(min.auc)) {
    if(!is.null(de[[1]]$AUC)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(AUC>min.auc))
      de <- de[unlist(lapply(de, nrow))>0]
      
    } else {
      warning("AUC column lacking in the DE results - recalculate with append.auc=TRUE")
    }
  }
  if(!is.null(min.specificity)) {
    if(!is.null(de[[1]]$Specificity)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(Specificity>min.specificity))
      de <- de[unlist(lapply(de, nrow))>0]
      
    } else {
      warning("Specificity column lacking in the DE results - recalculate append.specificity.metrics=TRUE")
    }
  }
  
  if(!is.null(min.precision)) {
    if(!is.null(de[[1]]$Precision)) {
      de <- lapply(de,function(x) x %>% dplyr::filter(Precision>min.precision))
      de <- de[unlist(lapply(de, nrow))>0]
      
    } else {
      warning("Precision column lacking in the DE results - recalculate append.specificity.metrics=TRUE")
    }
  }
  
  
  print('track')
  
  
  #de <- lapply(de,function(x) x%>%arrange(-Precision)%>%head(n.genes.per.cluster))
  if(n.genes.per.cluster==0) { # want to show only expliclty specified genes
    if(is.null(additional.genes)) stop("if n.genes.per.cluster is 0, additional.genes must be specified")
    additional.genes.only <- TRUE;
    n.genes.per.cluster <- 30; # leave some genes to establish cluster association for the additional genes
  } else {
    additional.genes.only <- FALSE;
  }
  
  de <- lapply(de,function(x) x%>%dplyr::arrange(!!rlang::parse_expr(ordering))%>%head(n.genes.per.cluster))
  de <- de[unlist(lapply(de, nrow))>0]
  
  gns <- lapply(de,function(x) as.character(x$Gene)) %>% unlist
  sn <- function(x) setNames(x,x)
  expl <- lapply(de,function(d) do.call(rbind,lapply(sn(as.character(d$Gene)),function(gene) conos:::getGeneExpression(con,gene))))
  
  # place additional genes
  if(!is.null(additional.genes)) {
    genes.to.add <- setdiff(additional.genes,unlist(lapply(expl,rownames)))
    x <- setdiff(genes.to.add,conos:::getGenes(con)); if(length(x)>0) warning('the following genes are not found in the dataset: ',paste(x,collapse=' '))
    
    age <- do.call(rbind,lapply(sn(genes.to.add),function(gene) conos:::getGeneExpression(con,gene)))
    # for each gene, measure average correlation with genes of each cluster
    acc <- do.call(rbind,lapply(expl,function(og) rowMeans(cor(t(age),t(og)),na.rm=T)))
    acc <- acc[,apply(acc,2,function(x) any(is.finite(x))),drop=F]
    acc.best <- na.omit(apply(acc,2,which.max))
    
    for(i in 1:length(acc.best)) {
      gn <- names(acc.best)[i];
      expl[[acc.best[i]]] <- rbind(expl[[acc.best[i]]],age[gn,,drop=F])
    }
    if(additional.genes.only) { # leave only genes that were explictly specified
      expl <- lapply(expl,function(d) d[rownames(d) %in% additional.genes,,drop=F])
      expl <- expl[unlist(lapply(expl,nrow))>0]
      
    }
  }
  
  
  exp <- do.call(rbind,expl)
  # limit to cells that were participating in the de
  exp <- na.omit(exp[,colnames(exp) %in% names(na.omit(groups))])
  
  if(order.clusters) {
    # group clusters based on expression similarity (of the genes shown)
    xc <- do.call(cbind,tapply(1:ncol(exp),groups[colnames(exp)],function(ii) rowMeans(exp[,ii,drop=F])))
    hc <- hclust(as.dist(2-cor(xc)),method='ward.D2')
    groups <- factor(groups,levels=hc$labels[hc$order])
    expl <- expl[levels(groups)]
    # re-create exp (could just reorder it)
    exp <- do.call(rbind,expl)
    exp <- na.omit(exp[,colnames(exp) %in% names(na.omit(groups))])
  }
  
  # transform expression values
  x <- t(apply(as.matrix(exp), 1, function(xp) {
    qs <- quantile(xp,c(1-expression.quantile,expression.quantile))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp-min(xp);
    xpr <- diff(range(xp));
    if(xpr>0) xp <- xp/xpr;
    xp
  }))
  
  
  o <- order(groups[colnames(x)])
  x=x[,o]
  
  
  annot <- data.frame(clusters=groups[colnames(x)],row.names = colnames(x))
  
  if(!is.null(column.metadata)) {
    if(is.data.frame(column.metadata)) { # data frame
      annot <- cbind(annot,column.metadata[colnames(x),])
    } else if(is.list(column.metadata)) { # a list of factors
      annot <- cbind(annot,data.frame(do.call(cbind.data.frame,lapply(column.metadata,'[',rownames(annot)))))
    } else {
      warning('column.metadata must be either a data.frame or a list of cell-named factors')
    }
  }
  annot <- annot[,rev(1:ncol(annot)),drop=FALSE]
  
  if(is.null(column.metadata.colors))  {
    column.metadata.colors <- list();
  } else {
    if(!is.list(column.metadata.colors)) stop("column.metadata.colors must be a list in a format accepted by HeatmapAnnotation col argument")
    # reorder pallete to match the ordering in groups
    if(!is.null(column.metadata.colors[['clusters']])) {
      if(!all(levels(groups) %in% names(column.metadata.colors[['clusters']]))) {
        stop("column.metadata.colors[['clusters']] must be a named vector of colors containing all levels of the specified cell groups")
      }
      column.metadata.colors[['clusters']] <- column.metadata.colors[['clusters']][levels(groups)]
    }
  }
  
  # make sure cluster colors are defined
  if(is.null(column.metadata.colors[['clusters']])) {
    uc <- unique(annot$clusters);
    column.metadata.colors$clusters <- setNames(rainbow(length(uc)),uc)
  }
  
  tt <- unlist(lapply(expl,nrow)); 
  rannot <- setNames(rep(names(tt),tt),unlist(lapply(expl,rownames)))
  #names(rannot) <- rownames(x);
  rannot <- rannot[!duplicated(names(rannot))]
  rannot <- rannot[names(rannot) %in% rownames(x)]
  rannot <- data.frame(clusters=factor(rannot,levels=names(expl)))
  
  x=unique(x)
  
  dim(rannot)
  dim(annot)
  dim(x)
  
  if(remove.duplicates) { x <- x[!duplicated(rownames(x)),] }
  
  # draw heatmap
  ha <- ComplexHeatmap::HeatmapAnnotation(df=annot,border=border,col=column.metadata.colors,show_legend=show.cluster.legend)
  

  if(show.gene.clusters) { 
    ra <- ComplexHeatmap::HeatmapAnnotation(df=rannot,which='row',show_annotation_name=FALSE, show_legend=FALSE, border=border,col=column.metadata.colors)
  } else { ra <- NULL } 
  
  #ComplexHeatmap::Heatmap(x, col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_column_names=FALSE, top_annotation=ha , left_annotation=ra, column_split=groups[colnames(x)], row_split=rannot[,1], row_gap = unit(0, "mm"), column_gap = unit(0, "mm"), border=T,  ...);
  if(split) {
    ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size), column_split=groups[colnames(x)], row_split=rannot[,1], row_gap = unit(split.gap, "mm"), column_gap = unit(split.gap, "mm"), ...);
  } else {
    ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size), ...);
    #ha <- ComplexHeatmap::Heatmap(x, name='expression', col=pal, cluster_rows=FALSE, cluster_columns=FALSE, show_row_names=is.null(labeled.gene.subset), show_column_names=FALSE, top_annotation=ha , left_annotation=ra, border=border,  show_heatmap_legend=show_heatmap_legend, row_names_gp = grid::gpar(fontsize = row.label.font.size));
    
  }
  if(!is.null(labeled.gene.subset)) {
    if(is.numeric(labeled.gene.subset)) {
      # select top n genes to show
      labeled.gene.subset <- unique(unlist(lapply(de,function(x) x$Gene[1:min(labeled.gene.subset,nrow(x))])))
    }
    gene.subset <- which(rownames(x) %in% labeled.gene.subset)
    labels <- rownames(x)[gene.subset];
    ha <- ha + ComplexHeatmap::rowAnnotation(link = ComplexHeatmap::anno_mark(at = gene.subset, labels = labels, labels_gp = grid::gpar(fontsize = row.label.font.size)))
    
  }
  
  if(return.details) {
    return(list(ha=ha,x=x,annot=annot,rannot=rannot,expl=expl,pal=pal,labeled.gene.subset=labeled.gene.subset))
  } else {
    return(ha)
  }
  
}


