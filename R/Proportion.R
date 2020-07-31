samToFrac=function(anoSample,fraction){
  inter=intersect(names(anoSample),names(fraction))
  tab=unique(data.frame('sample'=anoSample[inter],'frac'=fraction[inter]))
  rownames(tab)=tab$sample

  tab2=tab[,2]
  names(tab2)=rownames(tab)

  return(tab2)
}


proportion.box=function(anoCell,anoSample,fraction,appname,cname){

  ano2=data.frame('Cell'=anoCell[cname],'SampleType'=anoSample[cname])

  # Annotation vs sample
  tmp2 <- acast(ano2, Cell ~ SampleType, fun.aggregate=length)
  head(tmp2)
  # Normalise for the number of cells in each library
  tmp3 <- (sweep(tmp2, 2, colSums(tmp2), FUN='/'))
  tmp4 <- melt(tmp3)
  head(tmp4)
  names(tmp4) <- c('annot', 'sample','pc.of.sample')
  head(tmp4)

  p=ggplot(tmp4, aes(x=annot, fill=sample, y = pc.of.sample)) + geom_bar(stat='identity', position='fill') + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(appname,'.porpotion.png',sep=''),plot=p)


  tab=Toch(samToFrac(anoSample,fraction))

  tmp4$dtype=tab[as.character(tmp4$sample)]

  df=tmp4
  p <- ggplot(na.omit(df),aes(x=annot,y=pc.of.sample,dodge=dtype,fill=dtype))+geom_boxplot(notch=FALSE,outlier.shape=NA)  +  geom_point(position = position_jitterdodge(jitter.width=0.1),color=adjustcolor(1,alpha=0.3),aes(pch=dtype))+
    theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(angle = 90, hjust = 0.5))  +xlab("") +ylab("fraction of total cells")


  png(file=paste(appname,'.fraction.png',sep=''),width=600,height=400)
  print(p)
  dev.off();

  return(df)


}








