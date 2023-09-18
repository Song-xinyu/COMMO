timestart<-Sys.time()
#返回聚类结果，哪一列属于哪一类
ClusterPack=function(Data,cluster)
{
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Input:Cluster result, whose rank is the same as column names of "Data".
  #Output:Cluster output ordered by clustering assign(1,1...,1,2,2,...,2,3,3,...,3,...).
  temp=matrix(data = cluster,ncol=1)
  rownames(temp)=colnames(Data)
  temp=as.matrix(sort(temp[,1],na.last=T))
  return(temp)
}
#导入R包
library(magrittr)
library(dplyr)
#读入数据，行列转换
data<-read.table("D:/R/Rscripts/test20000×330.txt",row.names = 1, header = TRUE)
tdata<-t(data)
#设置聚类数
cluster_number=100
#用kmeans聚类
cluster_result=kmeans(data, cluster_number, nstart = 25)
#聚类结果
cluster_output=ClusterPack(Data=tdata,cluster=cluster_result$cluster)
#连接数据和聚类结果成新表格
join_data = merge(data, cluster_output, by = 0)
#按照类排序
join_data_order=join_data[order(join_data$V1),] 

class <- unique(join_data_order[,'V1'])
#导出总体聚类csv文件
result=join_data_order[,c(1,ncol(join_data_order))]
write_csv(result, 'D:/R/Rscripts/result/total.csv',append = TRUE)

#按照聚类结果分别写入txt文件
for (i in class){
  file_name = paste('module',i,'.txt',sep = '')
  file_path = paste('D:/R/Rscripts/result/',file_name,sep = '')
  file_bizcircle <- join_data_order %>% filter(V1 == i)
  aa=file_bizcircle[1:ncol(file_bizcircle)-1]
  write.table(aa,file_path,sep="\t" , row.names = FALSE ,quote = FALSE)
  
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime) 


