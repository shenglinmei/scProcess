

GOanalysis=function(markers,n){

  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]


  for(function_type in c("BP", "CC", "MF")){

    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)



    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)


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

}




GOanalysis.term=function(markers,n){

  ENTREZID=unlist(mget(markers, org.Hs.egSYMBOL2EG, ifnotfound=NA))
  ENTREZID=ENTREZID[!is.na(ENTREZID)]

  allr=list()

  for(function_type in c("BP", "CC", "MF")){

    param <- new("GOHyperGParams", geneIds=ENTREZID,
                 #universe=universe,
                 annotation="org.Hs.eg.db", ontology=function_type,pvalueCutoff=0.1,
                 conditional=FALSE, testDirection="over")
    hyp <- hyperGTest(param)
    sumTable <- summary(hyp)



    david=sumTable[1:20,]
    david$Pvalue=-log(david[,2])
    termNumber=nrow(david)

    david$genes=apply(david,1,function(x) { paste(names(ENTREZID[ENTREZID %in% get(as.character(x[1]),org.Hs.egGO2ALLEGS)]),collapse = ',' ) } )

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


