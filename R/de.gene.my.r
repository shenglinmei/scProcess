


getPerCellTypeDE2=function (con.obj, groups = NULL, sample.groups = NULL, cooks.cutoff = FALSE, 
                            ref.level = NULL, min.cell.count = 10, independent.filtering = FALSE, 
                            n.cores = 1, cluster.sep.chr = "<!!>", return.details = TRUE) 
{
  conos:::validatePerCellTypeParams(con.obj, groups, sample.groups, 
                                    ref.level, cluster.sep.chr)
  dat <- conos:::rawMatricesWithCommonGenes(con.obj, sample.groups)
  aggr2=lapply(dat,collapseCellsByType2, groups = groups, min.cell.count = min.cell.count) %>% 
    rbindDEMatrices2(cluster.sep.chr = cluster.sep.chr)
  
  gc()
  de.res <- conos:::papply(sn(levels(groups)), function(l) {
    tryCatch({
      cm <- aggr2[, conos:::strpart(colnames(aggr2), cluster.sep.chr, 
                                    2, fixed = TRUE) == l]
      meta <- data.frame(sample.id = colnames(cm), group = as.factor(unlist(lapply(colnames(cm), 
                                                                                   function(y) {
                                                                                     y <- conos:::strpart(y, cluster.sep.chr, 1, fixed = TRUE)
                                                                                     names(sample.groups)[unlist(lapply(sample.groups, 
                                                                                                                        function(x) any(x %in% y)))]
                                                                                   }))))
      if (!ref.level %in% levels(meta$group)) 
        stop("The reference level is absent in this comparison")
      meta$group <- relevel(meta$group, ref = ref.level)
      if (length(unique(as.character(meta$group))) < 2) 
        stop("The cluster is not present in both conditions")
      
      dds1 <- DESeq2::DESeqDataSetFromMatrix(cm, meta, 
                                             design = ~group)
      dds1 <- DESeq2::DESeq(dds1)
      res1 <- DESeq2::results(dds1, cooksCutoff = cooks.cutoff, 
                              independentFiltering = independent.filtering)
      res1 <- as.data.frame(res1)
      res1 <- res1[order(res1$padj, decreasing = FALSE), 
                   ]
      if (return.details) {
        list(res = res1, cm = cm, sample.groups = sample.groups)
      }
      else {
        res1
      }
    }, error = function(err) NA)
  }, n.cores = n.cores)
  de.res
}


#de.info <- getPerCellTypeDE2(con, groups=as.factor(new.annot), sample.groups = samplegroups, ref.level='bm', n.cores=4)


rbindDEMatrices2 <- function(mats, cluster.sep.chr) {
  #mats=aggr2
  index=unlist(lapply(mats,function(x) nrow(x)))>0
  mats=mats[names(index[index])]
  mats <- lapply(names(mats), function(n) {
    rownames(mats[[n]]) <- paste0(n, cluster.sep.chr, rownames(mats[[n]]));
    return(as.matrix(mats[[n]]))
  })
  #tmp=t(do.call(rbind, mats))
  return(t(do.call(rbind, mats)))
}





library(dplyr)

sn=function(x) { names(x) <- x; return(x); }


collapseCellsByType2=function (cm, groups, min.cell.count) 
{
  #cm=t[[1]]
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







library(xtable)
library(htmltools)

is.error2 = function (x)  {
  !('res' %in% names(x))
}






strpart <- function (x, split, n, fixed = FALSE) {
  sapply(strsplit(as.character(x), split, fixed = fixed), "[",n)
}

saveDEasJSON2 <- function(de.results = NULL, saveprefix = NULL, gene.metadata = NULL) {
  ## ### DEVEL
  ## de.results <- all.percl.TvsW
  ## saveprefix <- 'json/'
  ## rm(de.results, saveprefix)
  ## ##
  ## Check input
  
  ## Find de instances that didn't work (usually because cell type is absent from one or more sample types)
  n.error <- sum(unlist(lapply(de.results, is.error2)))
  print(n.error)
  if(n.error > 0)
    cat("Warning: ", n.error,' of ', length(de.results) ,' results have returned an error; ignoring...\n')
  ## get the de results that worked
  de.results <- de.results[!unlist(lapply(de.results, is.error2))]
  print(names(de.results))
  
  
  if(is.null(de.results)) stop('de.results have not been specified')
  if(is.null(saveprefix)) stop('saveprefix has not been specified')
  
  if(is.null(gene.metadata)) {
    gene.metadata <- data.frame()
    all.genes <- unique(unlist(lapply(de.results, function(x) {
      if(!is.null(x)){
        rownames(as.data.frame(x$res))
      } else {
        NULL
      }
    })))
    gene.metadata <- data.frame(geneid=all.genes)
  } else {
    if(is.null(gene.metadata$gene.id)) stop("gene.metadata must contain $gene.id field")  
  }
  
  
  
  ## Generate structure and save JSON
  lapply(names(de.results), function(ncc) {
    print(ncc)
    res.celltype <- de.results[[ncc]]
    ## Get results table as df
    res.table <- as.data.frame(res.celltype$res)
    ## append gene names
    res.table$gene <- rownames(res.table)
    ## append singificance
    res.table$significant <- res.table$padj < 0.05
    res.table$log2FoldChange[is.na(res.table$log2FoldChange)] <- 0 
    ## Append Z scores and rowid
    res.table$Z <- qnorm(1 - (res.table$pval/2))
    res.table$Z[is.na(res.table$Z)] <- 0
    res.table$Za <- qnorm(1 - (res.table$padj/2))
    res.table$Za[is.na(res.table$Za)] <- 0
    res.table$Z <- res.table$Z  * sign(res.table$log2FoldChange)
    res.table$Za <- res.table$Za  * sign(res.table$log2FoldChange)
    res.table$rowid <- 1:nrow(res.table)
    ## match order to metadata table
    
    #mo <- match(as.character(gene.metadata$geneid),as.character(res.table$gene))
    ## drop gene id column
    keep.cols <- colnames(gene.metadata)[colnames(gene.metadata) != 'geneid']
    names(keep.cols) <- keep.cols
    res.table <- cbind(res.table, gene.metadata[as.character(res.table$gene),keep.cols,drop=FALSE])
    ## get names of all the genes
    all.genes <- rownames(res.table)
    ## Get the count matrix
    cm <- res.celltype$cm
    ## remove the cell type suffix
    ## TODO make separator a parameter
    colnames(cm) <- strpart(colnames(cm),'<!!>',1,fixed=TRUE)
    ## ilev entry (submatrices of cps)
    ilev <- lapply(res.celltype$sample.groups, function(sg) {
      ## In certain cases columns may be missing,skip
      sg <- sg[sg %in% colnames(cm)]
      print(sg)
      ## keep only cols of interest
      cm.tmp <- cm[,sg]
      ## convert to matrix
      cm.tmp <- as.matrix(cm.tmp)
      rownames(cm.tmp) <- rownames(cm)
      ## calculate cpm
      cpm <- sweep(cm.tmp, 2, apply(cm.tmp,2, sum), FUN='/')
      cpm <- log10(cpm * 1e6 + 1)
      ## 
      snames1 <- colnames(cpm)
      ## Put genes in order
      cpm <- cpm[all.genes,]
      colnames(cpm) <- NULL;
      rownames(cpm) <- NULL;
      ## return
      list(snames=snames1, val=as.matrix(cpm))
    })
    ## snames entry (samplenames)
    snames <- names(res.celltype$sample.groups)
    ## convert to json
    tojson <- list(
      res = res.table,
      genes = all.genes,
      ilev = ilev,
      snames = snames
    )
    y <- jsonlite::toJSON(tojson)
    ## File to save to 
    file <- paste0(saveprefix,make.names(ncc),'.json')
    ## create the json file
    write(y,file)
    NULL
  })
  invisible(NULL)
}




ExpressionFraction=function(dat,l,cm,groups){
  cname=names(groups[groups==l])
  dat_tmp=lapply(dat[conos:::strpart(colnames(cm), cluster.sep.chr, 1, fixed = TRUE)],function(x) as.matrix(x[intersect(rownames(x),cname),]))
  
  tmp2=do.call(rbind,dat_tmp)
  tmp2[tmp2>0]=1
  ratio=Matrix::colSums(tmp2)/nrow(tmp2)
  return(ratio)
}


