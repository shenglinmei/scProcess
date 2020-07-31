
DEheatmap3=function(allg,exp,conosCluster,appname,removeGene=NULL,num=20){
  
  allg=lapply(allg,function(x) x[1:num])
  
  allg=lapply(allg,function(x) x[!is.na(x)])
  
  
  rowano=as.data.frame(stack(allg)[,1:2])
  
  uq=unique(rowano[,1])
  rowano=rowano[match(uq,rowano[,1]),]
  rownames(rowano)=rowano[,1]
  
  rowano=data.frame('group'=rowano[,2],row.names = rowano[,1])
  head(rowano)
  
  
  markers=unique(do.call(c,allg))
  markers=intersect(markers,colnames(exp))
  conosCluster=conosCluster[intersect(names(conosCluster),rownames(exp))]
  
  
  x <- as.matrix(exp[names(conosCluster),markers])
  ## trim outliers
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.01,0.98))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
  x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  
  
  #ss=apply(x,2,function(x) sum(x))
  #x=x[,abs(ss)>0.1]
  
  
  o <- order(conosCluster[colnames(x)])
  o2 <- order(rowano[rownames(x),1])
  
  
  x=x[o2,o]
  
  annot <- data.frame(CellType=conosCluster[colnames(x)],row.names = colnames(x))
  
  pal <- colorRampPalette(c('navy','white','firebrick3'))(50)
  ## draw heatmap
  
  fout=paste(appname,'.DiffG.png')
  rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
  heat=pheatmap(x,annotation_row = rowano,cluster_cols=FALSE,annotation_col = annot,show_colnames = F,annotation_legend = TRUE, #,gaps_col =  1,
                cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =5,width=8,height=8.8,   #3.5*0.02*length(markers),
                breaks = c(seq(min(x),-0.01,length.out = 50),0.01,seq(0.1,2,length.out = 48),max(x)))
   return(list('x'=x,'rowano'=rowano,'colano'=annot))
  
}








DEcaculate2.term.new=function(lallg,exp,appname,conosCluster,removeGene=NULL,cutoff=3,num=100,GO=NULL){
  allg=NULL  
  
  lallg=lapply(lallg,function(x) x[1:num])
  
  lallg=lapply(lallg,function(x) x[!is.na(x)])
  
  
  
  
  folderN=paste(appname,'.DE',sep='')
  pwd=getwd()
  pwd2=paste(pwd,'/',folderN,'/',sep='')
  
  system(paste('mkdir ',folderN))
  
  setwd(pwd2)
  

  for(n in names(lallg)){
    markers=lallg[[n]]
    allg=cbind(allg,markers)
    markers=intersect(markers,colnames(exp))
    x <- as.matrix(exp[names(conosCluster),markers])
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
             cluster_rows = FALSE,color=rgb.palette(100),filename=fout,fontsize_row =4,width=7,height=8,   #3.5*0.02*length(markers),
             breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
    
    
    iid= paste(appname,'_cluster_',n,sep='')
    if(!is.null(GO)){
      
      figGO=GOanalysis.term(markers,iid) 

    }
  }
  colnames(allg)=names(lallg)
  
  
  write.table(allg,paste(appname,'_Diff_gene.xls',sep=''),sep='\t',col.names=T,row.names=F,quote=F)
  setwd(pwd)
}  




getDifferentialGenes2=function (cm,type = "counts", clusterType = NULL, groups = NULL,
          name = "customClustering", z.threshold = 3, upregulated.only = FALSE,
          verbose = FALSE)
{
  "Determine differentially expressed genes, comparing each group against all others using Wilcoxon rank sum test\n\n       - type data type (currently only default 'counts' is supported)\n\n       - clusterType optional cluster type to use as a group-defining factor\n\n       - groups explicit cell group specification - a named cell factor (use NA in the factor to exclude cells from the comparison)\n\n       - name name slot to store the results in\n\n       - z.threshold minimal absolute Z score (adjusted) to report\n\n       - upregulated.only whether to report only genes that are expressed significantly higher in each group\n\n       - verbose verbose flag\n\n       return a list, with each element of the list corresponding to a cell group in the provided/used factor (i.e. factor levels); Each element of a list is a data frame listing the differentially epxressed genes (row names), with the following columns: Z - adjusted Z score, with positive values indicating higher expression in a given group compare to the rest; M - log2 fold change; highest- a boolean flag indicating whether the expression of a given gene in a given vcell group was on average higher than in every other cell group; fe - fraction of cells in a given group having non-zero expression level of a given gene;\n\n       Examples:\n\n\n         result <- r$getDifferentialGenes(groups=r$clusters$PCA$multilevel);\n\n         str(r$diffgenes)"
  if (is.null(groups)) {
    if (is.null(clusterType)) {
      cols <- clusters[[type]][[1]]
    }
    else {
      cols <- clusters[[type]][[clusterType]]
      if (is.null(cols)) {
        stop("clustering ", clusterType, " for type ",
             type, " doesn't exist")
      }
    }
  }
  else {
    cols <- groups
  }

  if (!all(rownames(cm) %in% names(cols))) {
    warning("cluster vector doesn't specify groups for all of the cells, dropping missing cells from comparison")
  }
  valid.cells <- rownames(cm) %in% names(cols)[!is.na(cols)]
  if (!all(valid.cells)) {
    cm <- cm[valid.cells, ]
  }
  cols <- as.factor(cols[match(rownames(cm), names(cols))])
  cols <- as.factor(cols)
  if (verbose) {
    cat("running differential expression with ", length(levels(cols)),
        " clusters ... ")
  }
  lower.lpv.limit <- -100
  xr <- pagoda2:::sparse_matrix_column_ranks(cm)
  grs <- pagoda2:::colSumByFac(xr, as.integer(cols))[-1, , drop = F]
  xr@x <- numeric(length(xr@x)) + 1
  gnzz <- pagoda2:::colSumByFac(xr, as.integer(cols))[-1, , drop = F]
  group.size <- as.numeric(tapply(cols, cols, length))[1:nrow(gnzz)]
  group.size[is.na(group.size)] <- 0
  gnz <- (group.size - gnzz)
  zero.ranks <- (nrow(xr) - diff(xr@p) + 1)/2
  ustat <- t((t(gnz) * zero.ranks)) + grs - group.size * (group.size +
                                                            1)/2
  n1n2 <- group.size * (nrow(cm) - group.size)
  usigma <- sqrt(n1n2 * (nrow(cm) + 1)/12)
  usigma <- sqrt((nrow(cm) + 1 - (gnz^3 - gnz)/(nrow(cm) *
                                                  (nrow(cm) - 1))) * n1n2/12)
  x <- t((ustat - n1n2/2)/usigma)
  if (verbose) {
    cat("adjusting p-values ... ")
  }
  x <- matrix(qnorm(pagoda2:::bh.adjust(pnorm(as.numeric(abs(x)), lower.tail = FALSE,
                                    log.p = TRUE), log = TRUE), lower.tail = FALSE, log.p = TRUE),
              ncol = ncol(x)) * sign(x)
  rownames(x) <- colnames(cm)
  colnames(x) <- levels(cols)[1:ncol(x)]
  if (verbose) {
    cat("done.\n")
  }
  log.gene.av <- log2(Matrix::colMeans(cm))
  group.gene.av <- pagoda2:::colSumByFac(cm, as.integer(cols))[-1, ,
                                                     drop = F]/(group.size + 1)
  log2.fold.change <- log2(t(group.gene.av)) - log.gene.av
  f.expressing <- t(gnzz/group.size)
  max.group <- max.col(log2.fold.change)
  if (upregulated.only) {
    ds <- lapply(1:ncol(x), function(i) {
      z <- x[, i]
      vi <- which(z >= z.threshold)
      r <- data.frame(Z = z[vi], M = log2.fold.change[vi,
                                                      i], highest = max.group[vi] == i, fe = f.expressing[vi,
                                                                                                          i])
      rownames(r) <- rownames(x)[vi]
      r <- r[order(r$Z, decreasing = T), ]
      r
    })
  }
  else {
    ds <- lapply(1:ncol(x), function(i) {
      z <- x[, i]
      vi <- which(abs(z) >= z.threshold)
      r <- data.frame(Z = z[vi], M = log2.fold.change[vi,
                                                      i], highest = max.group[vi] == i, fe = f.expressing[vi,
                                                                                                          i])
      rownames(r) <- rownames(x)[vi]
      r <- r[order(r$Z, decreasing = T), ]
      r
    })
  }
  names(ds) <- colnames(x)

  return(ds)
}


#matCSC <- as(t(tmp), "dgCMatrix")
#f=getDifferentialGenes2(matCSC,groups=conosCluster)






GOanalysis.term.mm=function(markers,n){

  #markers=rownames(table.down)[1:100]

  library(org.Mm.eg.db)
  ENTREZID=unlist(mget(markers, org.Mm.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]

  allr=list()

  for(function_type in c("BP", "CC", "MF")){

    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Mm.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)



    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)

    david$genes=apply(david,1,function(x) { paste(names(ENTREZID[ENTREZID %in% get(as.character(x[1]),org.Mm.egGO2ALLEGS)]),collapse = ',' ) } )

    allr[[function_type]]=david
    write.table(david,paste(n,'_' ,function_type, ".xls", sep=""),sep='\t',col.names=T,row.names=F,quote=F)


    library(ggplot2)
    p1 <- ggplot(data=david, aes(x=Pvalue, y=Term, size=Count, colour = factor(david$Count)))
    p1 <- p1 + geom_point()
    p1 <- p1 + guides(color = FALSE)
    p1 <- p1 + theme(panel.grid.major = element_line(colour='blue'),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank())
    p1 <- p1 + xlab(paste("-log10(Pvalue)", sep="")) + ylab("")
    p1 <- p1 + labs(title=paste("DAVID:", function_type, sep=""))
    p1 <- p1 + theme(axis.text.x=element_text(size=10, face="plain", colour ='black'))
    p1 <- p1 + theme(axis.text.y=element_text(size=6, face="plain", colour ='black'))
    #p1=p1+theme(axis.title.y=element_text(size=10,face="plain",colour ='black'))
    #p1=p1+theme(legend.position="bottom")
    p1 <- p1 + xlim(min(david$Pvalue), max(david$Pvalue))
    #p1=p1+ylim(-10,+15)
    print(p1)
    ggsave(file=paste(n,'_' ,function_type, ".png", sep=""), scale=0.8, dpi=600, width = 7, height=1+0.25*termNumber)
  }

  saveRDS(allr,paste(n,'.GOterm.rds',sep=''))
  return(allr)
}





 collapseCellsByType2=function (cm, groups, min.cell.count)
  {
    #cm=dat2[[1]]
    #groups=neu4
    g1 <- groups[intersect(names(groups), rownames(cm))]
    t1 <- as.numeric(table(g1))
    names(t1) <- levels(g1)
    droplevels <- names(t1)[t1 < min.cell.count]
    g1.n <- names(g1)
    g1 <- as.character(g1)
    names(g1) <- g1.n
    g1[g1 %in% droplevels] <- NA
    g1 <- as.factor(g1)
    aggr <- Matrix.utils::aggregate.Matrix(cm[names(g1),], g1)
    aggr <- aggr[rownames(aggr) != "NA", ]
    return(aggr)
  }




getClusterCountMatrices2=function (samples,clustering = NULL, groups = NULL, common.genes = TRUE,
            omit.na.cells = TRUE)
  {
    "Estimate per-cluster molecule count matrix by summing up the molecules of each gene for all of the cells in each cluster.\n\n\n       Params:\n\n       - clustering: the name of the clustering that should be used\n       - groups: explicitly provided cell grouping\n\n       - common.genes: bring individual sample matrices to a common gene list\n       - omit.na.cells: if set to FALSE, the resulting matrices will include a first column named 'NA' that will report total molecule counts for all of the cells that were not covered by the provided factor.\n       \n\n       Return: a list of per-sample uniform dense matrices with rows being genes, and columns being clusters\n      "
    if (is.null(groups)) {
      groups <- getClusteringGroups(clusters, clustering)
    }
    groups <- as.factor(groups)
    matl <- lapply(samples, function(s) {
      m <- t(s)
      cl <- factor(groups[match(rownames(m), names(groups))],
                   levels = levels(groups))
      tc <- conos:::colSumByFactor(m, cl)
      if (omit.na.cells) {
        tc <- tc[-1, , drop = F]
      }
      t(tc)
    })
     if (common.genes) {
      gs <- unique(unlist(lapply(matl, rownames)))
      matl <- lapply(matl, function(m) {
        nm <- matrix(0, nrow = length(gs), ncol = ncol(m))
        colnames(nm) <- colnames(m)
        rownames(nm) <- gs
        mi <- match(rownames(m), gs)
        nm[mi, ] <- m
        nm
      })
    }
    matl
}


