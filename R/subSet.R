
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
