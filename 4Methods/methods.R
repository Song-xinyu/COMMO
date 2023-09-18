####Flame-cluster method####
Flame_method=function(Data,cluster_result,ClusterNum=NULL, FLAMENum=NULL)
{
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))#Normalization for feature
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result$cluster)
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)
  return(result)
}

####K-means method####
Kmeans_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  #  getwd()
  library(magrittr)
  library(dplyr)
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))#Normalization for feature
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  print(cluster_number)
  cluster_result=kmeans(t(Data), cluster_number, nstart = 25) #Using the function of kmeans to calculate the cluster
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result$cluster) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)
  return(result)
}

####Hclust cluster method####
Hclust_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  library(magrittr)
  library(dplyr)
  library(readr)
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))#Normalization for feature
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  print(cluster_number)
  d=dist(t(Data), method = "euclidean")
  hc=hclust(d, method = "complete" )
  #hc=hclust(as.dist(1-cor(Data)),"ave") #Using the function of kmeans to calculate the cluster 修改点--------
  cluster_result=cutree(hc,k=cluster_number)
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  #survival_result=SurvivalPack(Survival_Information=Survival_Information,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)#删除survival_result=survival_result-------
  return(result)
}


####Spectral cluster method####
SpectralClust_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  library(SNFtool)
  library(magrittr)
  library(dplyr)
  library(readr)
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))#Normalization for feature
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  Distaa = (dist2(as.matrix(t(Data)),as.matrix(t(Data))))^(1/2)#修改点--------
  similarity=affinityMatrix(Distaa, 20, 0.5)#修改点--------
  cluster_result=spectralClustering(affinity=similarity,K=cluster_number) #Using the function of kmeans to calculate the cluster
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  #survival_result=SurvivalPack(Survival_Information=Survival_Information,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)#删除survival_result=survival_result----------------
  return(result)
}

####Agglomerative cluster method####
Agglomerative_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  library(cluster)
  library(plyr)
  require(plyr)
  library(dplyr)
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))#Normalization for feature
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  hc=agnes(t(Data), method = "complete")
  cluster_result=cutree(hc,k=cluster_number)
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)#删除survival_result=survival_result----------------
  return(result)
}


####ICA cluster method####
ICA_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  library(fastICA)
  library(readr)
  library(fdrtool)
  Data_rep=Data#Back up for signature illustration
  datExpr=t(Data)
  Data=t(as.matrix(scale(t(Data))))
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  
  a=fastICA(t(Data), cluster_number, alg.typ = "parallel", fun = "logcosh", alpha = 1,method = "R", row.norm = F, maxit = 1000, tol = 0.0001, verbose = TRUE)
  genes=NULL
  results=matrix(ncol=2)
  #colnames(results) = c("col1","col2")
  #功能点：每次循环需要添加上一轮循环的结果：gene cluster_result
  for (i in 1:ncol(a$S)){
    data = a$S[,i]
    fdr = fdrtool(data,statistic='normal',cutoff.method="fndr",plot=FALSE)   
    gene = names(fdr$pval[fdr$pval<0.05])#提取每列(每个样本)pval<0.05的基因-----不安全，如果一个没有会报错
    #join_data = datExpr[which(rownames(datExpr)%in% gene),]
    genes = append(genes,gene)
    result = matrix(nrow=length(gene),ncol=2)
    result[,1] = gene
    result[,2] = i
    results = rbind(results,result)
  }
  results = results[!duplicated(results[,1]),]#去重的基因
  rownames(results) = results[,1]#重新命名列
  results=as.matrix(results[-1,-1])#去除空值
  clusters = as.integer(results[,1])#取出聚类数
  cluster_output=matrix(data=clusters, ncol=1)#新建聚类输出
  rownames(cluster_output)=rownames(results)#匹配基因名
  genes=unique(genes)
  join_data = datExpr[which(rownames(datExpr)%in% genes),]#获得样本基因信息
  heatmap_result=SampleHeatmapPack(Data=t(join_data), ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)
  return(result)
}

####NMF method####
NMF_method_v2=function(Data, ClusterNum=NULL, FLAMENum=NULL)
{
  library(NMF)
  library(magrittr)
  library(dplyr)
  library(readr)
  Data_rep=Data#Back up for signature illustration
  Data=t(as.matrix(scale(t(Data))))
  Data_rank=apply(Data, 2, rank)
  Data_rank=scale(Data_rank, scale = TRUE, center = FALSE) #The number of NMF matrix must be greater than 0, so we can not use the parameter of center.
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  print(cluster_number)
  Res=nmf(Data_rank, cluster_number)
  cluster_result=NULL
  for (i in 1:ncol(coef(Res)))
  {
    cluster_result=c(cluster_result, which.max(coef(Res)[,i]))
  }
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)
  return(result)
}

####SOM method####
SOM_method_v2=function(Data, Method="kmeans", ClusterNum=NULL, FLAMENum=NULL)
{
  library(kohonen)
  library(magrittr)
  library(dplyr)
  Data_rep=Data#Back up for signature illustration
  Data=as.matrix(scale(t(Data)))
  #construct SOM and cluster samples into 5*5 subtypes
  som_grid=somgrid(xdim=10, ydim=10, topo = "hexagonal")
  resultsom=supersom(Data, grid=som_grid, keep.data = T)
  som_model_code_class=data.frame(name=rownames(Data), code_class=resultsom$unit.classif)
  som_data=as.matrix(as.data.frame(resultsom$codes))#codebook vector for 5*5 SOM point
  if(is.null(FLAMENum))# FLANME传参条件
  { 
    cluster_number=ClusterNum #Determine the optimal cluster number
  }
  else
  {
    cluster_number=FLAMENum
  }
  print(cluster_number)
  cluster_result=kmeans(som_data, cluster_number, nstart = 25) #Using the function of kmeans to calculate the cluster
  som_model_code_class=cbind(som_model_code_class,som_model_code_class[,2])
  colnames(som_model_code_class)[3]="cluster"
  som_cluster=cluster_result$cluster
  # 
  # if (0)
  # {
  #   if (is.null(FLAMENum))# FLANME传参条件
  #   { 
  #     cluster_number=FLAMENum #Determine the optimal cluster number
  #   }
  #   else
  #   {
  #     cluster_number=ClusterNum
  #   }
  #   hc=hclust(as.dist(1-cor(t(som_data))),"ave") #Using the function of kmeans to calculate the cluster
  #   cluster_result=cutree(hc,k=cluster_number)
  #   som_model_code_class=cbind(som_model_code_class,som_model_code_class[,2])
  #   colnames(som_model_code_class)[3]="cluster"
  #   som_cluster=cluster_result
  # }
  
  for (i in 1:length(unique(som_cluster)))
  {
    som_model_code_class[which(som_model_code_class[,2] %in% which(som_cluster==i)),3]=i
  }
  Data=t(Data)
  cluster_result=som_model_code_class[,3]
  cluster_output=ClusterPack(Data=Data,cluster=cluster_result) 
  heatmap_result=SampleHeatmapPack(Data=Data, ClusterPackResult=cluster_output)
  network_result=SampleNetworkPack(Similarity=heatmap_result,ClusterPackResult=cluster_output)
  #survival_result=SurvivalPack(Survival_Information=Survival_Information,ClusterPackResult=cluster_output)
  signature_result=SignaturePack(Data=Data_rep,ClusterPackResult=cluster_output)
  result=list(cluster_number=cluster_number,cluster_output=cluster_output,heatmap_result=heatmap_result,
              network_result=network_result,signature_result=signature_result)#删除survival_result=survival_result
  return(result)
}
