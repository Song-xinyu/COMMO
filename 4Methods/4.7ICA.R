timestart<-Sys.time()
#导入数据
datExpr <- read.table("D:/R/Rscripts/test1500×200.txt", row.names = 1,header = TRUE,encoding ='utf-8')
library(fastICA)

#fsatICS聚类，记得改聚类数
a <- fastICA(datExpr,15, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = F, maxit = 1000,
             tol = 0.0001, verbose = TRUE)

library(readr)
library(fdrtool)

#提取聚类后每个样本的pvalue<0.05的基因，分别放在新txt里，并导出总体csv表格
for (i in 1:ncol(a$S)){
  data <- a$S[,i]
  fdr <- fdrtool(data,statistic='normal',cutoff.method="fndr",plot=FALSE)   
  gene <- names(fdr$pval[fdr$pval<0.05])#提取每列(每个样本)pval<0.05的基因
  join_data<-datExpr[which(rownames(datExpr)%in% gene),]
  file_name = paste('module',i,'.txt',sep = '')
  file_path = paste('D:/R/Rscripts/result/',file_name,sep = '')
  
  write.table(join_data,file_path,sep="\t" , row.names = T ,quote = FALSE)
  
  result <- matrix(nrow=length(gene),ncol=2)
  result[,1] <- gene
  result[,2] <- i
  result <- as.data.frame(result)
  
  write_csv(result, 'D:/R/Rscripts/result/total.csv',append = TRUE)
  
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 