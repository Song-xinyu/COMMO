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
library(kohonen)
library(magrittr)
library(dplyr)
#读入数据，行列转换
data<-read.table("D:/R/Rscripts/test1500×200.txt",row.names = 1, header = TRUE)
tdata<-t(data)
Datas=as.matrix(scale(data))

som_grid=somgrid(xdim=10, ydim=10, topo = "hexagonal")
#SOM聚类
resultsom=supersom(Datas, grid=som_grid, keep.data = T)
som_model_code_class=data.frame(name=rownames(Datas), code_class=resultsom$unit.classif)
som_data=as.matrix(as.data.frame(resultsom$codes))

############hc和Kmean都可以运行
# cluster_number=50
# hc=hclust(as.dist(1-cor(t(som_data))),"ave") #Using the function of hc to calculate the cluster
# cluster_result=cutree(hc,k=cluster_number)
# som_model_code_class=cbind(som_model_code_class,som_model_code_class[,2])
# colnames(som_model_code_class)[3]="cluster"
# som_cluster=cluster_result
###################
cluster_number=15
cluster_result=kmeans(som_data, cluster_number, nstart = 25) #Using the function of kmeans to calculate the cluster
som_model_code_class=cbind(som_model_code_class,som_model_code_class[,2])
colnames(som_model_code_class)[3]="cluster"
som_cluster=cluster_result$cluster

###############

for (i in 1:length(unique(som_cluster)))
{
  som_model_code_class[which(som_model_code_class[,2] %in% which(som_cluster==i)),3]=i
}

Data1=t(Datas)###########??????
cluster_result=som_model_code_class[,3]

cluster_output=ClusterPack(Data=Data1,cluster=cluster_result)
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



