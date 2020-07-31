#







#gs=c('CD3E','CD8A','CD8B','CD4','KLRD1','FCGR3A','FGFBP2','GZMA','GZMB','GZMK','CCR7','SELL','LEF1','CCR7','CXCR3','CCR6')

#gs=c('SOX11','SOX4','CDKN3','CDK1','MKI67','STMN2','TH','PHOX2B','NPY')

geneList_exp=function(con,gs,exp,fout,alpha =0.6,size=0.4){
  print(setdiff(gs,rownames(exp)))
  gs=intersect(gs,rownames(exp))
  lcol=4
  lrow=celling(length(gs)/4)
  lis=list()
  for (gene in gs){
    print(gene)
    
    a=con$plotGraph(colors =exp[gene,nname],alpha =alpha,plot.na=F,size=size,title=gene)
    lis[[gene]]=a
  }
  b=  cowplot::plot_grid(plotlist=lis, ncol=lrow, nrow=lcol)
  ggsave(fout,b,width = 2*lrow,height=2.05*lcol)
}





piechart=function(sample){
  # pie chart
  ratio=(table(sample)/length(sample))*100
  
  # Create Data
  data <- data.frame(
    cell=paste(names(ratio),': ',round(ratio,2),'%',sep=''),
    value=as.numeric(ratio)
  )
  
  # Basic piechart
  a1=ggplot(data, aes(x="", y=value, fill=cell)) +
    geom_bar(stat="identity", width=1) +
    coord_polar("y", start=0)+
    scale_fill_manual(values=as.character(annot.palf(length(unique(sample.palf)))))+
    theme_void()
  return(a1)
}




qplot <- function(g, con.obj, ann) {
  #cat(g,' ')
  x <- lapply(con.obj$samples[5:10],function(r) { if(g%in% colnames(r$counts)) { r$counts[,g] } else { return(NULL) } })
  if(length(unlist(x))<1) stop('gene ',g,' is not found')
  df <- data.frame(val=unlist(x),cell=unlist(lapply(x,names)))
  df$cluster <- ann[match(df$cell,names(ann))]
  df <- na.omit(df)
  
  mv <- max(tapply(df$val,df$cluster,quantile,p=0.8),tapply(df$val,df$cluster,mean))*1.5
  p <- ggplot(df,aes(x=cluster,y=val,color=cluster))+geom_boxplot(outlier.shape = NA)+ stat_summary(fun.data=mean_se,geom="pointrange", color="black")+ylab(g)+ggtitle(g)+guides(colour=FALSE)+ theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+coord_cartesian(ylim=c(0,mv));
  p
}



