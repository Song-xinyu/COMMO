timestart<-Sys.time()
library(cluster)
library(plyr)
require(plyr)
library(dplyr)
#导入数据
datExpr <- read.table("D:/R/Rscripts/test20000×330.txt", row.names = 1,header = TRUE,encoding ='utf-8')
data <- scale(datExpr)#标准化
#用agnes聚类
hc <- agnes(data, method = "complete")
sub_grp <- cutree(hc, k = 100)#改聚类数！！！
table(sub_grp)
#
join_data = mutate(datExpr,cluster = sub_grp)
#按照类排序
join_data_order=join_data[order(join_data$cluster),] 
result=data.frame(join_data_order[,ncol(join_data)],row.names =rownames(join_data_order) )
#导出总体聚类csv文件
write.csv(result, 'D:/R/Rscripts/result/total.csv')

class <- unique(join_data_order[,'cluster'])
#按照聚类结果分别写入txt文件
for (i in class){
  file_name = paste('module',i,'.txt',sep = '')
  file_path = paste('D:/R/Rscripts/result/',file_name,sep = '')
  file_bizcircle <- join_data_order %>% filter(cluster == i)
  aa=file_bizcircle[1:ncol(file_bizcircle)-1]
  write.table(aa,file_path,sep="\t" , row.names = T ,quote = FALSE)
  
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 
