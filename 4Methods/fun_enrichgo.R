# GO富集分析
# gene:仅包含Entrez ID(txt/csv)
# pvalueCutoff/qvalueCutoff:通路显著性阈值
# minGSSize/maxGSSize:通路中包含的最小/最大基因数
# USER_DATA:GO富集分析背景文件,包含GeneID、GO_ID、GO_term和Category四列

fun_enrich_internal_go = function(gene,pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                                  minGSSize = 10,maxGSSize = 500,USER_DATA){
  library(dplyr)
  library(stringr)
  
  N = length(unique(USER_DATA$GeneID)) # 参与GO term的所有基因数目
  list.intersect = as.character(intersect(gene,unique(USER_DATA$GeneID))) # 差异基因与全部基因取交集
  n = length(list.intersect) # 差异基因与全部基因的交集基因数
  
  if(n > 0){
    tmp.df = USER_DATA
    tmp.df$TF = tmp.df$GeneID %in% list.intersect # 标记参与GO term的差异基因
    bp.stat = data.frame(table(tmp.df$GO_ID)) # 统计每个term的基因数
    colnames(bp.stat)[1] = "GO_ID"
    bp.stat = bp.stat %>% filter(Freq >= minGSSize & Freq <= maxGSSize) # 保留满足条件的term
    tmp.stat = subset(tmp.df,TF==TRUE)
    tmp.stat = data.frame(table(tmp.stat$GO_ID))
    colnames(tmp.stat)[1] = "GO_ID"
    colnames(tmp.stat)[2] = "Count"
    #tmp.stat = data.frame(tmp.df %>% group_by(GO_ID) %>% summarize(Count = sum(TF))) # 每个term富集到的差异基因数
    tmp.stat = tmp.stat %>% filter(Count > 0) # 保留富集到差异基因的term
    keep.pw = sort(intersect(bp.stat$GO_ID,tmp.stat$GO_ID)) # 保留差异基因富集到的和满足条件的term
    if(length(keep.pw) > 0){
      # 计算p值
      path.term = c()
      pvalue = c()
      for (i in keep.pw) {
        one.set = tmp.df[tmp.df$GO_ID %in% i,]
        M = length(one.set$GeneID) # 注释到某个term的基因数
        x = sum(list.intersect %in% one.set$GeneID) # 富集到某个term的差异基因数
        one.pvalue = phyper(x-1,M,N-M,n,lower.tail = FALSE) # 超几何检验计算p值
        path.term = append(path.term,i)
        pvalue = append(pvalue,one.pvalue)
      }
      phyper.res = data.frame(GO_ID = path.term,pvalue = pvalue)
      phyper.res = phyper.res %>% arrange(pvalue) # 按照p值升序排序
      phyper.res$p.adjust = p.adjust(phyper.res$pvalue,method = "fdr") # 计算fdr值
      
      # 计算q值
      qobj = tryCatch(qvalue::qvalue(p = phyper.res$pvalue,lambda = 0.05,pi0.method = "bootstrap"),error = function(e) NULL)
      if(class(qobj) == "qvalue"){qvalue = qobj$qvalues}else{qvalue = NA}
      phyper.res$qvalue = qvalue
      
      phyper.merge = merge(phyper.res,merge(bp.stat,tmp.stat,by = "GO_ID")) # 合并term的总/差异基因数
      colnames(phyper.merge)[5] = "Total"
      
      # 计算GeneRatio/BgRatio/Fold Enrichment等并添加到结果
      phyper.merge$GeneRatio = paste0(phyper.merge$Count,"/",n)
      phyper.merge$BgRatio = paste0(phyper.merge$Total,"/",N)
      # phyper.merge$FoldEnrichment = (phyper.merge$Count/n)/(phyper.merge$Total/N)
      
      # 添加通路名称
      df = subset(USER_DATA,subset = GO_ID %in% phyper.merge$GO_ID,select = -GeneID)
      df = distinct(df)
      finall = merge(df,phyper.merge)
      
      # 添加GO term包含的差异基因ID
      gene.set = tmp.df %>% filter(TF == TRUE) # 保留富集到差异基因的记录
      gene.set = subset(gene.set,subset = GO_ID %in% finall$GO_ID,select = c(GeneID,GO_ID)) # 保留分析结果的记录
      library(purrr)
      library(data.table)
      fun_paste0 = function(x){
        tmp = data.frame(GO_ID = x$GO_ID[1],geneID = paste0(x$GeneID,collapse = "/"))
        return(tmp)
      }
      gene.df = gene.set %>% split(.$GO_ID) %>% map(fun_paste0) %>% rbindlist()
      finall = merge(finall,gene.df,by = "GO_ID")
      finall = finall[order(finall$p.adjust),] # 按照p.adjust升序排序
      colnames(finall)[1:3] = c("ID","Description","ONTOLOGY")
      
      # 选取p.adjust和qvalue符合阈值的结果
      Over = finall[finall$pvalue <= pvalueCutoff,]
      Over = Over[Over$p.adjust <= pvalueCutoff,]
      if(!any(is.na(Over$qvalue))){Over = Over[Over$qvalue <= qvalueCutoff,]}
      return(Over)
    }else{
      return(Over = data.frame())
      message("--> fun_enrich_internal_go: No GO term can be matched....")
    }
  }else{
    return(Over = data.frame())
    message("--> fun_enrich_internal_go: No gene can be used....")
  }
}

fun_enrich_go = function(gene,pvalueCutoff = 0.05,qvalueCutoff = 0.2,ont = "ALL",
                         minGSSize = 10,maxGSSize = 500,USER_DATA){
  library(dplyr)
  library(stringr)
  
  # 获取背景文件
  USER_DATA = read.table(USER_DATA,header = T,sep = "\t",quote = "",stringsAsFactors = F)
  
  if(length(gene) > 0){
    go.bp = USER_DATA %>% filter(Category == "BP") %>%
      fun_enrich_internal_go(gene = gene,pvalueCutoff = pvalueCutoff,qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,maxGSSize = maxGSSize)
    go.cc = USER_DATA %>% filter(Category == "CC") %>%
      fun_enrich_internal_go(gene = gene,pvalueCutoff = pvalueCutoff,qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,maxGSSize = maxGSSize)
    go.mf = USER_DATA %>% filter(Category == "MF") %>%
      fun_enrich_internal_go(gene = gene,pvalueCutoff = pvalueCutoff,qvalueCutoff = qvalueCutoff,
                             minGSSize = minGSSize,maxGSSize = maxGSSize)
    go = rbind(go.bp,go.cc,go.mf)
    # 返回结果
    if(ont == "BP"){return(go.bp)}else if(ont == "CC"){return(go.cc)}else if(ont == "MF"){return(go.mf)}else if(ont == "ALL"){return(go)}
  }else{
    message("--> fun_enrich_go: No gene can be used....")
  }
}