

# ncdcon$getDatasetPerCell()

sn=function(x) { names(x) <- x; return(x); }


mergeMatrix=function(raw.mats){
  
  genelists <- lapply(raw.mats, function(x) rownames(x))
  str(genelists)
  commongenes <- Reduce(intersect,genelists)
  
  
  matrices2 <- mapply(function(m, name) {
    #   colnames(m) <- paste(name, colnames(m), sep='_');
    m[commongenes,]
  },
  raw.mats,
  names(raw.mats))
  cellNum=unlist(lapply(matrices2, function(x) ncol(x)))
  print('#commongenes')
  print(length(commongenes))
  
  print('#cell numbers')
  print(cellNum)
  bigM2 <- Reduce(cbind, matrices2)
  
  return(bigM2)
}



#  Run robust ranking
library(RobustRankAggreg)
robust_ranking=function(res,decreasing=NULL) {
  
  DataOrder = list()
  DataOrderName = list()
  nset=ncol(res)
  for (i in seq(nset)){
    order_tmp =order(as.numeric(as.character(res[,i])),decreasing=F)
    if (!is.null(decreasing)){ order_tmp =order(as.numeric(as.character(res[,i])),decreasing=TRUE) }
    tmp_name = rownames(res[order_tmp,])
    print(i)
    print(tmp_name[1:4])
    print(order_tmp[1:4])
    DataOrderName[[i]]=tmp_name
    DataOrder[[i]] = rank(as.numeric(as.character(res[,i])))}
  
  glist = 	DataOrderName
  r = RobustRankAggreg::rankMatrix(glist) 
  roubustRanking = aggregateRanks(rmat = r)
  return(roubustRanking)
}

Toch=function(conosCluster){
  cluster2=as.character(conosCluster)
  names(cluster2)=names(conosCluster)
  return(cluster2)
}



outputVelocity=function(appname,group4,emb,spliced,unspliced){
  #  output(i,bridge.ano[tnname],con$embedding[tnname,],spliced,unspliced)
  obs=data.frame('name'=names(group4),'clusters'=group4)
  cell=names(group4)
  
  splicedf=spliced[,cell]
  unsplicedf=unspliced[,cell]
  gs=rownames(splicedf)
  
  Matrix::writeMM(splicedf,paste(appname,'.splicedf.mtx',sep=''))
  Matrix::writeMM(unsplicedf,paste(appname,'.unsplicedf.mtx',sep=''))
  write.csv(obs,paste(appname,'.obs.csv',sep=''),row.names=F)
  write.csv(emb,paste(appname,'.emb.csv',sep=''),row.names=F,col.names=F)
  write.csv(gs,paste(appname,'.gs.csv',sep=''),row.names=F,col.names=T)
}


disF=function(con,typefc,samplef){
  cm <- scon$getClusterCountMatrices(groups=typefc)
  samples=scon$samples
  
  ctdm <- lapply(cm,function(x) {apply(x,2,function(x) log10((x/pmax(1,sum(x)))*1e3+1 ) )} )
  
  ctdm2 <- lapply(sn(colnames(ctdm[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(ctdm,function(x) x[,ct]))
    tcm
  })
  
  
  table(names(typefc)==names(samplef))
  
  cct <- table(typefc,samplef)
  dis <- lapply(sn(colnames(cm[[1]])),function(ct) {
    tcm <- do.call(rbind,lapply(cm,function(x) x[,ct]))
    tcm <- t(tcm/pmax(1,rowSums(tcm)))
    tcd <- pagoda2:::jsDist(tcm); dimnames(tcd) <- list(colnames(tcm),colnames(tcm));
    # calculate how many cells there are
    attr(tcd,'cc') <- cct[ct,colnames(tcm)]
    tcd
  })

  #x[upper.tri(x)] <- NA; diag(x) <- NA;
  #df2 <- na.omit(melt(x))
  
  return(list('CPM'=ctdm2,'dist'=dis))
}



getMarkers=function(){
  gg='INOS|IL12|FCGR1A|FCGR1B|FCGR1C|CD80|CXCR10|IL23|CXCL9|CXCL10|CXCL11|CD86|IL1A|IL1B|IL6|TNFa|MHCII|CCL5|IRF5|IRF1|CD40|IDO1|KYNU|CCR7'
  M1=strsplit(gg,split='|', fixed=TRUE)[[1]]
  
  gg='ARG1|ARG2|IL10|CD32|CD163|CD23|FCER2|CD200R1|PD-L2|PD-L1|MARCO|CSF1R|CD206|IL1RA|IL14R|CCL4|CCL13|CCL20|CCl17|CCL18|CCl22|CCL24|LYVE1|VEGFA|VEGFB|VEGFC|VEGFD|EGF|CTSA|CTSB|CTSC|CTSD|TGFB1|TGFB2|TGFB3|MMP14|MMP19|MMP9|CLEC7A|WNT7B|FASL|TNSF12|TNSF8CD276|VTCN1|MSR1|FN1|IRF4'
  M2=strsplit(gg,split='|', fixed=TRUE)[[1]]
  
  gs='ATF3,BCL2A1,CCL20,CXCL2,DUSP2,EREG,IL1B,PDE4B,PTGS2,PTX3,TNF,TNFAIP3,IL6,G0S2,S100A8,IL8,IL1RN,TREM1,OSM,NLRP3,IL10,AQP9'
  TIM=strsplit(gs,split=',', fixed=TRUE)[[1]]
  
  f1='/d0-mendel/home/meisl/bin/data/score/GSE5099_Macropage_VS_Mono0h_downMono.txt'
  f2='/d0-mendel/home/meisl/bin/data/score/GSE5099_Macropage_VS_Mono0h_upMacro.txt'
  f3='/d0-mendel/home/meisl/bin/data/score/GSE8286_Macropage_VS_Mono0h_downMono.txt'
  f4='/d0-mendel/home/meisl/bin/data/score/GSE8286_Macropage_VS_Mono0h_upMacro.txt'
  
  gs=read.csv(f1,sep='\t',header=F)
  GSE5099_Mono=as.character(gs[,1])
  
  gs=read.csv(f2,sep='\t',header=F)
  GSE5099_Macro=as.character(gs[,1])
  
  gs=read.csv(f3,sep='\t',header=F)
  GSE8286_Mono=as.character(gs[,1])
  
  gs=read.csv(f4,sep='\t',header=F)
  GSE8286_Macro=as.character(gs[,1])
  
  
  gs=read.csv('/d0-mendel/home/meisl/bin/data/180209_cell_cycle_regev.txt',sep='\t',header=F)
  cycle=as.character(gs[,1])
 
  cytotoxicity=c('GZMA','GZMB','GZMM','GZMK','GZMH','PRF1','CD8A','CD8B')
  
  
  gs=read.csv('/d0-mendel/home/meisl/bin/data/Treg.activity.txt',sep='\t',header=F)
  TregActivity=as.character(gs[,1])
  
  
  ehaust3=read.csv('/d0-mendel/home/meisl/bin/data/exhausted.gene.sig',sep='\t',header=F)
  ehaust3.zeming=as.character(ehaust3[,1])
  
  
  ehaust1=read.csv('/d0-mendel/home/meisl/Workplace/BMME/b.newAno/cell/Exhasut.marker2',sep='\t',header=F)
  ehaust.regv=as.character(ehaust1[,1])
  
  tmp=read.csv('/d0-mendel/home/meisl/bin/data/cell_cycle_Scanpy.txt',sep='\t',header=T)
  
  
  dat=list('M1'=M1,
           'M2'=M2,
           'TIM'=TIM,
           'GSE5099_Mono'=GSE5099_Mono,
           'GSE5099_Macro'=GSE5099_Macro,
           'GSE8286_Mono'=GSE8286_Mono,
           'GSE8286_Macro'=GSE8286_Macro,
           'cycle'=cycle,
           'cytotoxicity'=cytotoxicity,
           'ehaust3.zeming'=ehaust3.zeming,
           'ehaust.regv'=ehaust.regv,
           'TregActivity'=TregActivity,
            'G1.S'=gsub(' ','',as.character(tmp[,2])) %>% .[.!=''],
           'S'=gsub(' ','',as.character(tmp[,3])) %>% .[.!=''],
           'G2.S'=gsub(' ','',as.character(tmp[,4])) %>% .[.!=''],
           'M'=gsub(' ','',as.character(tmp[,5])) %>% .[.!=''],
           'M.G1'=gsub(' ','',as.character(tmp[,6])) %>% .[.!='']
           )
  return(dat)
}


getSampleNamePerCell=function (samples) {
  cl <- lapply(samples, getCellNames)
  return(rep(names(cl), sapply(cl, length)) %>% stats::setNames(unlist(cl)) %>% 
           as.factor())
}

cluster.expression.distances2=function(con,groups=NULL,dist='JS',n.cores=con$n.cores,min.cluster.size=1,min.samples=1,max.n.cells=Inf,aggr=median,return.details=FALSE,use.aggregated.matrices=FALSE,use.single.cell.comparisons=FALSE,n.cells=100,lib.size.scale=1e6) {
  require(abind)
  if(is.null(groups)) {
    if(is.null(con$clusters)) stop('no groups specified and no clusterings found') 
    groups <- as.factor(con$clusters[[1]]$groups)
  } else {
    groups <- as.factor(groups)
  }
  
  valid.dists <- c('JS','cor','PCA');
  if(!dist %in% valid.dists) stop(paste('only the following distance types are supported:',paste(valid.dists,collapse=', ')))
  
  # limit the factor only to the cells that appear within the dataset
  groups <- droplevels(groups[names(groups) %in% rownames(con$embedding)])  
  
  if(use.aggregated.matrices) {
    # in this case, we simply add all the counts from all the clusters into a common matrix and simply calculate distances on that
    tc <- conos:::rawMatricesWithCommonGenes(con) %>% 
      lapply(conos:::collapseCellsByType, groups=as.factor(groups), min.cell.count=0) %>%
      abind(along=3) %>%
      apply(c(1,2),sum,na.rm=T)
    
    if(dist=='JS') {
      tc <- t(tc/pmax(1,rowSums(tc)))
      tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
    } else { # correlation distance
      tc <- log10(t(tc/pmax(1,rowSums(tc)))*lib.size.scale+1)
      tcd <- 1-cor(tc)
    }
    diag(tcd) <- 0;
    if(return.details) {
      return(list(dc=tc,mdist=tcd))
    } else {
      return(tcd)
    }
  }
  
  
  # determine distance matrices for each sample
  dcl <- conos:::papply(con$samples,function(s) {
    m <- s$misc$rawCounts[rownames(s$counts), ] # $counts cells can be size filtered ... 
    cl <- factor(groups[match(rownames(m),names(groups))],levels=levels(groups));
    tt <- table(cl);
    
    empty <- tt<min.cluster.size;
    
    if(use.single.cell.comparisons) { # sample individual cells and compare
      
      tcd <- lapply(1:n.cells,function(i) {
        scn <- unlist(tapply(names(cl),cl,function(x) sample(x,1)))
        tc <- as.matrix(m[na.omit(as.character(scn)),,drop=F])
        rownames(tc) <- names(scn)[!is.na(scn)]
        tc <- tc[match(levels(cl),rownames(tc)),]
        rownames(tc) <- levels(cl)
        if(dist=='JS') {
          tc <- t(tc/pmax(1,rowSums(tc)))
          tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
        } else if(dist=='PCA') {
          tcn <- s$reductions$PCA[match(scn,rownames(s$reductions$PCA)), ,drop=FALSE]; rownames(tcn) <- rownames(tc);
          tcd <- 1-cor(t(tcn))
        } else { # correlation distance
          tc <- log10(t(tc/pmax(1,rowSums(tc)))*lib.size.scale+1)
          tcd <- 1-cor(tc)
        }
        tcd[empty,] <- tcd[,empty] <- NA;
        tcd
      }) %>% abind(along=3) %>% apply(c(1,2),median,na.rm=T)
    } else { # aggregated clusters
      if(any(tt>max.n.cells)) { # need subsampling
        scn <- unlist(tapply(names(cl),cl,function(x) sample(x,min(max.n.cells,length(x)))))
        cl[!(names(cl) %in% scn)] <- NA; tt <- table(cl);
      }
      tc <- conos:::colSumByFactor(m,cl);
      tc <- tc[-1,,drop=F]  # omit NA cells
      if(dist=='JS') {
        tc <- t(tc/pmax(1,rowSums(tc)))
        tcd <- pagoda2:::jsDist(tc); dimnames(tcd) <- list(colnames(tc),colnames(tc));
      } else if(dist=='PCA') {
        # determine medians
        cl <- na.omit(cl);
        tcn <- do.call(rbind,tapply(names(cl),cl,function(nn) apply(s$reductions$PCA[nn,,drop=FALSE],2,median)))
        # correct for missing levels
        tcn <- tcn[match(levels(cl),rownames(tcn)),]; rownames(tcn) <- levels(cl);
        tcd <- 1-cor(t(tcn))
        #tcne <- matrix(NA,nrow=length(levels(cl)),ncol=ncol(tcn)); rownames(tcne) <- levels(cl); colnames(tcne) <- colnames(tcn); ml <- match(rownames(tcn),rownames(tcne)); tcne[ml,] <- tcn;
        #tcd <- 1-cor(t(tcne))
      } else { # correlation distance
        tc <- log10(t(tc/pmax(1,rowSums(tc)))*lib.size.scale+1)
        tcd <- 1-cor(tc)
      }
    }
    
    tcd[empty,] <- tcd[,empty] <- NA;
    diag(tcd) <- 0;
    # calculate how many cells there are
    attr(tcd,'cc') <- table(cl)
    tcd
  },n.cores=n.cores,mc.preschedule=T)
  
  nc <- lapply(dcl,attr,'cc')
  dc <- abind(dcl,along=3)
  
  
  # summarize across samples
  n.valid.obs <- apply(dc,c(1,2),function(x) sum(!is.na(x)))
  mdist <- apply(dc,c(1,2),aggr,na.rm=T)
  
  mdist[n.valid.obs<min.samples] <- NA;
  if(return.details) {
    return(list(dc=dc,nc=nc,mdist=mdist))
  } else {
    return(mdist)
  }
}
