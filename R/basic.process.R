read10xMatrixm=function (path) {
  matrixFile <- paste(path, "/matrix.mtx.gz", sep = "")
  genesFile <- paste(path, "/features.tsv.gz", sep = "")
  barcodesFile <- paste(path, "/barcodes.tsv.gz", sep = "")
  if (!file.exists(matrixFile)) {
    stop("Matrix file does not exist")
  }
  if (!file.exists(genesFile)) {
    stop("Genes file does not exist")
  }
  if (!file.exists(barcodesFile)) {
    stop("Barcodes file does not exist")
  }
  x <- as(Matrix::readMM(gzfile(matrixFile)), "dgCMatrix")
  genes <- read.table(gzfile(genesFile))
  rownames(x) <- genes[, 2]
  barcodes <- read.table(gzfile(barcodesFile))
  colnames(x) <- barcodes[, 1]
  invisible(x)
}





rawToPagoda=function(datraw,appname,tSNE=TRUE){

  genelists <- lapply(datraw, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  length(commongenes)


  matrices_raw <- mapply(function(mm, name) {
    mm[commongenes,]
  },
  datraw,
  names(datraw))


  bigM2 <- Reduce(cbind, matrices_raw)


  p2=runPagoda(bigM2,appname,n.cores = 12,get.tsne=tSNE)

  f1=paste(appname,'_p2combined.rds',sep='')
  saveRDS(p2,f1)


  return(p2)
}





runConos=function(datlp2,appname,n.cores=8){
  print(names(datlp2))
  con <- Conos$new(datlp2,n.cores=n.cores)
  con$buildGraph()
  con$findCommunities()


  f1=paste(appname,'_conos_clustering.pdf',sep='')
  p1=con$plotGraph()
  ggsave(f1,p1)

  f1=paste(appname,'_conos.rds',sep='')
  saveRDS(con,f1)

  con$embedGraph(method = 'UMAP')

  p=con$plotGraph()
  ggsave(paste(appname,'.umap.png',sep=''),p,height=7,width=7)


  saveRDS(con$embedding,paste(appname,'.umap.rds',sep=''))


  con$embedGraph(method="UMAP", n.cores=30,min.dist=1e-10,n.neighbors=50)

  p2a <- con$plotGraph(alpha=0.1,size=0.1,plot.na=T)

  saveRDS(con$embedding,paste(appname,'.umap2.rds',sep=''))


  ggsave(paste(appname,'.umap409.png',sep=''),p2a,height=7,width=7)


  return(con)


}




#  run pagoda apps
runPagoda=function(cd,appname='none', n.cores = 1, batch = NULL, n.odgenes = 3000, nPcs = 100,
                   k = 30, perplexity = 50, log.scale = TRUE, trim = 10, keep.genes = NULL,
                   min.cells.per.gene = 0, get.largevis = TRUE, get.tsne = TRUE,
                   make.geneknn = TRUE) {

  rownames(cd) <- make.unique(rownames(cd))
  p2 <- Pagoda2$new(cd, n.cores = n.cores, batch = batch, keep.genes = keep.genes,
                    trim = trim, log.scale = log.scale, min.cells.per.gene = min.cells.per.gene)


  pv=paste(appname,'.adjustVariance.pdf',sep='')
  pdf(pv)
  p2$adjustVariance(plot = T, gam.k = 10)
  dev.off()

  p2$calculatePcaReduction(nPcs = nPcs, n.odgenes = n.odgenes,
                           maxit = 1000)
  p2$makeKnnGraph(k = k, type = "PCA", center = TRUE, weight.type = "none",
                  n.cores = n.cores, distance = "cosine")
  p2$getKnnClusters(method = igraph::infomap.community, type = "PCA",
                    name = "infomap")
  p2$getKnnClusters(method = igraph::multilevel.community,
                    type = "PCA", name = "multilevel")

  p2$getKnnClusters(method = igraph::walktrap.community, type = 'PCA', name='walktrap')

  p2$getDifferentialGenes(type='PCA',verbose=T,clusterType='multilevel')


  p2$getEmbedding(type = 'PCA', embeddingType = 'tSNE', perplexity = 50)
  #    M <- 30
  #    p2$getEmbedding(type = 'PCA', embeddingType = 'largeVis', M = M, perplexity = perplexity, gamma = 1 / M, alpha = 1 )


  pdf1=paste(appname,'.tsn.pdf',sep='')
  pdf(pdf1)
  p2$plotEmbedding(type='PCA',embeddingType='tSNE',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (tSNE)')
  dev.off()



  return(p2)
}







fsubset=function(raw,anoCell,input_cell,cutoff=30){

  cname=names(anoCell[anoCell %in% input_cell])

  raw2=lapply(raw,function(x) x[,intersect(cname,colnames(x))])

  num=unlist(lapply(raw2,function(x) ncol(x)))


  print(num)

  num=num[ num>cutoff ]
  raw2=raw2[names(num) ]

  print(unlist(lapply(raw2,function(x) ncol(x))))

  return(raw2)
}





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

