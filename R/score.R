
getScore=function(aexp,gs,group,anoSample,anoType,magnitude.normal=NULL,rscore=NULL){
  # exp colname is gene 
  
  gs=intersect(gs,rownames(aexp))
  cname=names(group)
  exp=as.matrix(aexp[gs,cname])
  
  exp=exp[rowMeans(exp)!=0,] # remove unexpressed gene 
  exp=exp[,colSums(exp)!=0] 
  cname=intersect(cname,colnames(exp))
  
  if (!is.null(magnitude.normal)){
    exp=apply(exp,2,function(x) x/max(x))
  }

  score=rowMeans(t(exp))
  score_scale=scale(score)
  names(score_scale)=names(score)
  
  
  dat=data.frame('score'=score_scale,'cell'=group[cname],'type'=anoType[cname],'sample'=anoSample[cname])
  dat$type2=apply(dat,1,function(x) paste(x['cell'],'|',x['sample'],'|',x['type'],sep=''))
  
  tmp=tapply(dat$score,dat$type2,mean)
  index=match(names(tmp),dat$type2)

  dat2=data.frame('name'=dat$type2[index],'score'=tmp,'type'=dat$type[index],'cell'=dat$cell[index],'sample'=dat$sample[index])
  
  if (!is.null(rscore)){
    return(score_scale)
  }else{
    return(dat2)
  }
  
}





ScoreMatrix=function(gs,aexp,group,anoSample,anoType,plot=NULL,file=NULL){
  
 # aexp=Aexp
 #group=typef[typef %in% c('Mono1')]
 # cname=names(group)
  #anoSample=samplef[cname]
  #anoType=fractionf[cname]
  
  
  gs=intersect(gs,rownames(aexp))
  cname=names(group)
  exp=as.matrix(aexp[gs,cname])
  
  ttype=paste(anoSample,'|',anoType,sep='')
    
  mat=apply(exp,1,function(x) {
    tapply(x,ttype,mean)
  })
  
  
  
  x=mat
  x <- apply(x, 2, function(xp) {
    qs <- quantile(xp,c(0.05,0.95))
    xp[xp<qs[1]] <- qs[1]
    xp[xp>qs[2]] <- qs[2]
    xp
  })
 # x <- x[,(apply(x,2,sd) != 0)]
  x <- t(scale(x))
  ano=lapply(sn(colnames(x)),function(x) strsplit(x,'[|]')[[1]])
  sample=unlist(lapply(ano,function(x) x[1]))
  faction=unlist(lapply(ano,function(x) x[2]))
  annot = data.frame(sample,faction,row.names = colnames(x))
  
  
  if (!is.null(plot)){
    o <- order(annot$faction)
    rgb.palette <- colorRampPalette(c("blue","white","red"), space = "rgb" )
    aa=pheatmap(x[,o],annotation_col = annot,show_rownames = T,show_colnames = T,width=4+0.1*ncol(x),
                cluster_rows = TRUE, cluster_cols=FALSE,color=rgb.palette(100),filename=file,height=3+0.1*nrow(x),
                breaks = c(seq(min(x),-0.01,length.out = 50),seq(0.01,max(x),length.out = 50)))
  }
  
  return('ano'=annot,'exp'=x)
    
}
  
  
  
  




