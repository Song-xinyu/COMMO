# KEGG富集分析
# gene(Entrez ID,txt/csv)
# pvalueCutoff/qvalueCutoff:通路显著性阈值
# minGSSize/maxGSSize:通路中包含的最小/最大基因数
# USER_DATA:KEGG富集分析背景文件,包含path.id、entrez.id和path.name三列

fun_enrich_internal = function(gene,pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                               minGSSize = 10,maxGSSize = 500,USER_DATA){
  library(dplyr)
  
  # 获取背景文件
  USER_DATA = read.table(USER_DATA,header = T,sep = "\t",stringsAsFactors = F)
  
  N = length(unique(USER_DATA$entrez.id)) # 参与KEGG通路的所有基因数目
  list.intersect = as.character(intersect(gene,USER_DATA$entrez.id)) # 差异基因与参与通路的所有基因取交集
  n = length(list.intersect) # 差异基因与参与通路的所有基因的交集基因数
  
  if(n > 0){
    tmp.df = USER_DATA
    tmp.df$TF = tmp.df$entrez.id %in% list.intersect # 标记参与通路的差异基因
    bp.stat = data.frame(table(tmp.df$path.id)) # 统计每条通路的基因数
    colnames(bp.stat)[1] = "path.id"
    bp.stat = bp.stat %>% filter(Freq >= minGSSize & Freq <= maxGSSize) # 保留满足条件的通路
    tmp.stat = subset(tmp.df,TF==TRUE)
    tmp.stat = data.frame(table(tmp.stat$path.id))
    colnames(tmp.stat)[1] = "path.id"
    colnames(tmp.stat)[2] = "Count"
    #tmp.stat = data.frame(tmp.df %>% group_by(path.id) %>% summarize(Count = sum(TF))) # 统计每条通路富集到的差异基因数
    tmp.stat = tmp.stat %>% filter(Count > 0) # 保留富集到差异基因的通路
    keep.pw = sort(intersect(bp.stat$path.id,tmp.stat$path.id)) # 保留满足条件的差异基因富集到的通路
    if(length(keep.pw) > 0){
      # 计算p值
      path.term = c()
      pvalue = c()
      for (i in keep.pw) {
        one.set = tmp.df[tmp.df$path.id %in% i,]
        M = length(one.set$entrez.id) # 注释到某条通路的基因数
        x = sum(list.intersect %in% one.set$entrez.id) # 富集到某条通路的差异基因数
        one.pvalue = phyper(x-1,M,N-M,n,lower.tail = FALSE) # 超几何检验计算p值
        path.term = append(path.term,i)
        pvalue = append(pvalue,one.pvalue)
      }
      phyper.res = data.frame(path.id = path.term,pvalue = pvalue)
      phyper.res = phyper.res %>% arrange(pvalue) # 按照p值升序排序
      phyper.res$p.adjust = p.adjust(phyper.res$pvalue,method = "fdr") # 计算fdr值
      
      # 计算q值
      qobj = tryCatch(qvalue::qvalue(p = phyper.res$pvalue,lambda = 0.05,pi0.emthod = "bootstrap"),error = function(e) NULL)
      if(class(qobj) == "qvalue"){qvalue = qobj$qvalues}else{qvalue = NA}
      phyper.res$qvalue = qvalue
      
      # 合并通路中的总/差异基因数目
      phyper.merge = merge(phyper.res,merge(bp.stat,tmp.stat,by = "path.id"))
      colnames(phyper.merge)[5] = "Total"
      
      # 计算GeneRatio/BgRatio/FoldEnrichment等并添加到结果
      phyper.merge$GeneRatio = paste0(phyper.merge$Count,"/",n)
      phyper.merge$BgRatio = paste0(phyper.merge$Total,"/",N)
      # phyper.merge$FoldEnrichment = (phyper.merge$Count/n)/(phyper.merge$Total/N)
      
      # 添加通路名称
      df = subset(USER_DATA,subset = path.id %in% phyper.merge$path.id,select = -entrez.id)
      df = distinct(df)
      finall = merge(df,phyper.merge)
      
      # 添加KEGG通路包含的差异基因ID
      gene.set = tmp.df %>% filter(TF == TRUE) # 保留富集到差异基因的记录
      gene.set = subset(gene.set,subset = path.id %in% finall$path.id,select = c(path.id,entrez.id))
      library(purrr)
      library(data.table)
      fun_paste0 = function(x){
        tmp = data.frame(path.id = x$path.id[1],geneID = paste0(x$entrez.id,collapse = "/"))
        return(tmp)
      }
      gene.df = gene.set %>% split(.$path.id) %>% map(fun_paste0) %>% rbindlist()
      finall = merge(finall,gene.df)
      finall = finall[order(finall$p.adjust),] # 按照p.adjust升序排序
      colnames(finall)[1:2] = c("ID","Description")
      
      # 选取p.adjust和qvalue符合阈值的结果
      Over = finall[finall$pvalue <= pvalueCutoff,]
      Over = Over[Over$p.adjust <= pvalueCutoff,]
      if(!any(is.na(Over$qvalue))){Over = Over[Over$qvalue <= qvalueCutoff,]}
      return(Over)
    }else{
      return(Over = data.frame())
      message("--> fun_enrich_internal: No KEGG pathway can be matched....")
    }
  }else{
    return(Over = data.frame())
    message("--> fun_enrich_internal: No gene can be used....")
  }
}