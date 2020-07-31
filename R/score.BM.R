# @param  aexp expression matrix; row is cell, coloum is gene
# @param sample.cell ratio of sampling cell size,  range from 0 to 1 
# @param  sample.gs  ratio of sampling gene number ,  range from 0 to 1 
Signature_score=function(anoCell,m2,aexp,fraction,anoSample,min.num.cell=10,magnitude.normal = NULL,sample.cell=NULL,sample.gs=NULL){
  m2=intersect(m2,colnames(aexp))
  cname=names(anoCell)
  
  if (!is.null(sample.cell)){
    sample.cell=ceiling(length(cname)*sample.cell)
    print(sample.cell)
    cname=sample(cname,sample.cell)
  }
  
  if (!is.null(sample.gs)){
    sample.gs=ceiling(length(gs)*sample.gs)
    m2=sample(m2,sample.gs)
  }
  
  
  exp=as.matrix(aexp[cname,m2])
  
  mmax=apply(exp,2,function(x) max(x))
  exp=exp[,mmax>0]
  
  if (!is.null(magnitude.normal)) {
    exp = apply(exp, 2, function(x) x/max(x))
  }  
  m2score=rowMeans(as.matrix(exp))
  summary(m2score)
  
  m2score_scale=scale(m2score)
  names(m2score_scale)=names(m2score)
  print(summary(m2score))
  
  cname=names(m2score)
  
  
  ttype=paste(anoCell[cname],'|',anoSample[cname],'|',fraction[cname],sep='')
  names(ttype)=cname
  
  if (!is.null(min.num.cell)){
    sel= table(ttype) %>% .[.>min.num.cell] %>% names()
    ttype=ttype[ttype %in% sel]
  }
 # print(table(ttype))
  
  m2score=m2score[names(ttype)]
  score.sample=tapply(m2score,ttype,mean)
  
  ano=lapply(sn(names(score.sample)),function(x) strsplit(x,'[|]')[[1]])
  cell=unlist(lapply(ano,function(x) x[1]))
  anoSample=unlist(lapply(ano,function(x) x[2]))
  faction=unlist(lapply(ano,function(x) x[3]))
  
  dat2=data.frame('score'=score.sample,'cell'=cell,'fraction'=faction,'sample'=anoSample)
  dat2$name=paste(as.character(dat2[,'cell']),as.character(dat2[,'sample']))
  dat2$fraction=ordered(as.factor(dat2$fraction), levels = levels(fraction))   
  dat2$cell=ordered(as.factor(dat2$cell), levels = levels(anoCell))   
  
#  dat2=dat2[dat2$fraction=='Tumor',]
  return(dat2)
}
  




drawBoxplot=function(tmp,name2,myeloid.col,limHeight=1.5,dsize=3,sigl=NULL,height=3.3,width=2.5,ysize=11,xsize=11){
  # limHeight=1.45
  t=tapply(tmp$score , tmp$cell,median)
  nn=names(t)[which.max(t)]
  
  # filter significant pairs 
  sig=compare_means(score ~ cell,  data = tmp)
  write.table(sig,paste(name2,'.pvalue.xls',sep=''),col.names=T,row.names=F,quote=F,sep='\t')
  sig=sig[sig$p.signif!='ns',]

  if (!is.null(sigl)){
    
    sig=sig[sig$group1 %in% sigl,]
    sig=sig[sig$group2 %in% sigl,]
  }
  
  siglis=split(sig, seq(nrow(sig)))
  pair=lapply(siglis,function(x) as.character(x[,2:3]))
  
  p1 <- ggplot(tmp, aes(x=cell,fill=cell,y=score)) + geom_boxplot(outlier.shape = -1,width=0.5,position=position_dodge(width=0.1)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")  + ylab(name2)
  p1=p1+ geom_point(data = tmp,,color=adjustcolor(1,alpha=0.3), position = position_jitterdodge(0.3)) +
    stat_compare_means(comparisons = pair,label = "p.signif",hide.ns=TRUE,size=dsize,tip.lengt=0.01,bracket.size =0.3,label.y.npc = "bottom")  #+ # Add pairwise comparisons p-value
  
  p1=p1+ylim(c(min(tmp$score)*.7,max(tmp$score)*limHeight))
  #p1=p1+ geom_point(data = tmp,color=adjustcolor(1,alpha=0.3),fill='grey', size = 1, shape = 21)
  p1=p1+ theme(legend.position="none")
  p1=p1+theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p1=p1+theme( axis.text.y = element_text(angle = 45, hjust = 0.5,color = "black"),axis.text.x=element_text(size=xsize,color = "black"),axis.title.y = element_text(size = ysize,color = "black"))
 
  if (!is.null(myeloid.col)){
    p1=p1+scale_fill_manual(values=myeloid.col)
  }
  
  
  ggsave(paste(name2,'.score.pvalue.pdf',sep=''),p1,w=width,h=height)
  saveRDS(tmp,paste(name2,'.dat.rds',sep=''))
  
  #p1 <- ggplot(tmp, aes(x=cell,fill=cell,y=score)) + geom_boxplot(outlier.shape = -1,width=0.5,position=position_dodge(width=0.1)) +theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + xlab("")  + ylab(name2)
  #p1=p1+ geom_point(data = tmp,,color=adjustcolor(1,alpha=0.3), position = position_jitterdodge(0.3)) + scale_fill_manual(values=myeloid.col)
  #p1=p1+ theme(legend.position="none")+ylim(c(min(tmp$score)*.7,max(tmp$score)*1.4))
  #p1
  #ggsave(paste(name2,'.score.pdf',sep=''),p1,w=width,h=height)
  #print(p1) 
  return(p1)  
}




getExp_Sample=function(aexp,gene,group,anoSample,fraction,min.num.cell=5,scale=NULL){
  # exp colname is gene

  cname=intersect(names(group),names(fraction))

  gexp=aexp[gene,cname]

  ttype=paste(group[cname],'|',anoSample[cname],'|',fraction[cname],sep='')
  names(ttype)=cname

  if (!is.null(min.num.cell)){
    sel= table(ttype) %>% .[.>min.num.cell] %>% names()
    ttype=ttype[ttype %in% sel]
  }

  cname=names(ttype)
  Mexp=tapply(gexp[cname],ttype,mean)

  if (!is.null(scale)){
    score_Mexp = scale(Mexp)
    names(score_Mexp) = names(Mexp)
    Mexp=score_Mexp
  }
  ano=lapply(sn(names(Mexp)),function(x) strsplit(x,'[|]')[[1]])
  group=unlist(lapply(ano,function(x) x[1]))
  anoSample=unlist(lapply(ano,function(x) x[2]))
  fraction=unlist(lapply(ano,function(x) x[3]))

  dat2=data.frame('score'=Mexp,'group'=group,'sample'=anoSample,'fraction'=fraction)
  return(dat2)
}

