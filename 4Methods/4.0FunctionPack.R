####A pack of functions about analysis####

####1 Select optimal cluster number####
####1.1 Select optimal cluster number using 5 index (for kmeans and hclust)####
OptClustNum=function(Data,Method="kmeans")
{
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Input:Method,"kmeans"or"hc".
  #Output:Optimal cluster number.
  library(NbClust)
  seed=8715
  index=c("ch","db","silhouette","dunn")
  results=NULL
  for(i in 1:length(index))
  {
    if (Method=="kmeans")
    {
      temp=try(NbClust(t(Data), min.nc=2, max.nc=8, method="kmeans", index=index[i]),silent = TRUE)
    }
    if (Method=="hc")
    {
      temp=try(NbClust(t(Data), min.nc=2, max.nc=8, method="average", index=index[i]),silent = TRUE)
    }
    if('try-error' %in% class(temp))
    {
      next
    }else{
      results=c(results,temp$Best.nc[1])
    }
  }
  temp=as.matrix(table(results))
  optimal_cluster_num=max(as.numeric(rownames(temp)[which.max(temp)]))
  return(optimal_cluster_num)
}
####1.2 Select optimal cluster number using CDF (for ikmeans and ihclust)####
OptClustNum_Fori=function(ConsensusClusterPlus_result)
{
  #Input:Result from ConsensusClusterPlus package.
  #Output:Optimal cluster number.
  Kvec=2:8
  x1=0.1
  x2=0.9 # threshold defining the intermediate sub-interval
  PAC=rep(NA,length(Kvec))
  names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK
  for(i in Kvec){
    M = ConsensusClusterPlus_result[[i]]$consensusMatrix
    Fn = ecdf(M[lower.tri(M)])
    PAC[i-1] = Fn(x2) - Fn(x1)
  }#end for i
  optimal_cluster_num = Kvec[which.min(PAC)]
  return(optimal_cluster_num)
}
####1.3 Select optimal cluster number using gap (for spectral clustering)####
OptClustNum_ForSpec=function(Similarity)
{
  #Input:Sample similarity matrix.
  #Output:Optimal cluster number.
  library(SNFtool)
  return(estimateNumberOfClustersGivenGraph(Similarity, 2:8)[[1]])
}
####1.4 Select optimal cluster number using cophenetic correlation coefficient (for NMF)####
OptClustNum_ForNMF=function(Data)
{
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Output:Optimal cluster number.
  library(NMF)
  #Res_multi=nmf(Data, 2:8, nrun=10, seed=8715) 
  Res_multi=nmf(Data, 2:8)
  coph=Res_multi$measures$cophenetic
  coph_diff=NULL
  for (i in 2:length(coph)) 
  {
    coph_diff=c(coph_diff, coph[i-1]-coph[i])
  }
  optimal_cluster_num=which.max(coph_diff)+1
  return(optimal_cluster_num)
}

####2 Pack cluster result for visualization####
####2.1 Pack cluster result####
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
####2.2 Pack CMS cluster result####
SampleCMSPack_old=function(SparseCrossHeatmap,ClusterPackResult,ClustMethod)
{
  SparseCrossHeatmap=SparseCrossHeatmap[rownames(ClusterPackResult),]
  sample_cluster_matrix=matrix(data = 0,nrow = length(unique(ClusterPackResult)),ncol = ncol(SparseCrossHeatmap))
  colnames(sample_cluster_matrix)=colnames(SparseCrossHeatmap)
  rownames(sample_cluster_matrix)=unique(ClusterPackResult)
  for (i in 1:ncol(sample_cluster_matrix))
  {
    temp=ClusterPackResult[which(SparseCrossHeatmap[,i]==1)]
    for (j in 1:nrow(sample_cluster_matrix))
    {
      if (is.na(rownames(sample_cluster_matrix)[j]))
      {
        sample_cluster_matrix[j,i]=length(which(is.na(temp)))
      }else
      {
        sample_cluster_matrix[j,i]=sum(rownames(sample_cluster_matrix)[j]==temp,na.rm = T)
      }
    }
    rm(temp)
  }
  p_value_mapping=matrix(data = NA,nrow = length(ClustMethod)+1,ncol = 2)
  colnames(p_value_mapping)=c("frequency","p_value")
  p_value_mapping[,1]=0:length(ClustMethod)
  for (i in 1:nrow(p_value_mapping))
  {
    N=length(ClustMethod)*nrow(sample_cluster_matrix)
    n=length(ClustMethod)
    M=length(ClustMethod)
    k=p_value_mapping[i,1]
    p_value_mapping[i,2]=phyper(k-1,M, N-M, n, lower.tail=FALSE)
  }
  sample_cluster_matrix_pvalue=sample_cluster_matrix
  for (i in 1:nrow(p_value_mapping))
  {
    sample_cluster_matrix_pvalue[which(sample_cluster_matrix==p_value_mapping[i,1])]=p_value_mapping[i,2]
  }
  PAVLUE_CUTTOF=0.05
  index1=which(sample_cluster_matrix_pvalue<=PAVLUE_CUTTOF)
  index2=which(sample_cluster_matrix_pvalue>PAVLUE_CUTTOF)
  sample_cluster_matrix_pvalue[index1]=1
  sample_cluster_matrix_pvalue[index2]=0
  result=matrix(data = NA,nrow = ncol(sample_cluster_matrix_pvalue),ncol = 1)
  rownames(result)=colnames(sample_cluster_matrix_pvalue)
  for (i in 1:nrow(sample_cluster_matrix_pvalue))
  {
    index=which(sample_cluster_matrix_pvalue[i,]==1)
    if (length(index)>0)
    {
      result[index,1]=i
    }
  }
  result=as.matrix(sort(result[,1],na.last=T))
  return(result)
}

SampleCMSPack=function(SparseCrossHeatmap,ClusterPackResult,ClustMethod)
{
  SparseCrossHeatmap=SparseCrossHeatmap[rownames(ClusterPackResult),]
  sample_cluster_matrix=matrix(data = 0,nrow = length(unique(ClusterPackResult)),ncol = ncol(SparseCrossHeatmap))
  colnames(sample_cluster_matrix)=colnames(SparseCrossHeatmap)
  rownames(sample_cluster_matrix)=unique(ClusterPackResult)
  for (i in 1:ncol(sample_cluster_matrix))
  {
    temp=ClusterPackResult[which(SparseCrossHeatmap[,i]==1)]
    for (j in 1:nrow(sample_cluster_matrix))
    {
      if (is.na(rownames(sample_cluster_matrix)[j]))
      {
        sample_cluster_matrix[j,i]=length(which(is.na(temp)))
      }else
      {
        sample_cluster_matrix[j,i]=sum(rownames(sample_cluster_matrix)[j]==temp,na.rm = T)
      }
    }
    rm(temp)
  }
  sample_cluster_matrix_pvalue=sample_cluster_matrix
  for (i in 1:ncol(sample_cluster_matrix))
  {
    for (j in 1:nrow(sample_cluster_matrix))
    {
      if (!is.na(rownames(sample_cluster_matrix)[j]))
      {
        N=nrow(ClusterPackResult)
        n=length(ClustMethod)
        M=length(which(ClusterPackResult==rownames(sample_cluster_matrix)[j]))
        k=sample_cluster_matrix[j,i]
        sample_cluster_matrix_pvalue[j,i]=phyper(k-1,M, N-M, n, lower.tail=FALSE)
      }else
      {
        sample_cluster_matrix_pvalue[j,i]=1
      }
    }
  }
  PVALUE_CUTOFF=0.05
  for (i in 1:ncol(sample_cluster_matrix_pvalue))
  {
    index1=which(sample_cluster_matrix_pvalue[,i]>=PVALUE_CUTOFF)
    if (length(index1)>0)
    {
      sample_cluster_matrix_pvalue[index1,i]=1
    }
    if (min(sample_cluster_matrix_pvalue[,i])==1)
    {
      sample_cluster_matrix_pvalue[,i]=0
    }else{
      index=which.min(sample_cluster_matrix_pvalue[,i])
      sample_cluster_matrix_pvalue[index,i]=1
      sample_cluster_matrix_pvalue[-index,i]=0
      rm(index)
    }
  }

  result=matrix(data = NA,nrow = ncol(sample_cluster_matrix_pvalue),ncol = 1)
  rownames(result)=colnames(sample_cluster_matrix_pvalue)
  for (i in 1:nrow(sample_cluster_matrix_pvalue))
  {
    index=which(sample_cluster_matrix_pvalue[i,]==1)
    if (length(index)>0)
    {
      result[index,1]=i
    }
  }
  result=as.matrix(sort(result[,1],na.last=T))
  return(result)
}

SampleCMSPackFortwothree=function(SparseCrossHeatmap,ClusterPackResult,ClustMethod)
{
  SparseCrossHeatmap=SparseCrossHeatmap[rownames(ClusterPackResult),]
  sample_cluster_matrix=matrix(data = 0,nrow = length(unique(ClusterPackResult)),ncol = ncol(SparseCrossHeatmap))
  colnames(sample_cluster_matrix)=colnames(SparseCrossHeatmap)
  rownames(sample_cluster_matrix)=unique(ClusterPackResult)
  for (i in 1:ncol(sample_cluster_matrix))
  {
    temp=ClusterPackResult[which(SparseCrossHeatmap[,i]==1)]
    for (j in 1:nrow(sample_cluster_matrix))
    {
      if (is.na(rownames(sample_cluster_matrix)[j]))
      {
        sample_cluster_matrix[j,i]=length(which(is.na(temp)))
      }else
      {
        sample_cluster_matrix[j,i]=sum(rownames(sample_cluster_matrix)[j]==temp,na.rm = T)
      }
    }
    rm(temp)
  }
  index=which(sample_cluster_matrix>0.5*length(ClustMethod))
  sample_cluster_matrix[index]=1
  sample_cluster_matrix[-index]=0
  rm(index)
  
  result=matrix(data = NA,nrow = ncol(sample_cluster_matrix),ncol = 1)
  rownames(result)=colnames(sample_cluster_matrix)
  for (i in 1:nrow(sample_cluster_matrix))
  {
    index=which(sample_cluster_matrix[i,]==1)
    if (!is.na(rownames(sample_cluster_matrix)[i]))
    {
      if (length(index)>0)
      {
        result[index,1]=i
      }
    }else
    {
      result[index,1]=NA
    }
  }
  result=as.matrix(sort(result[,1],na.last=T))
  return(result)
}

####3 Pack sample similarity heatmap for visualization####
####3.1 Pack sample similarity heatmap for feature*sample matrix####
SampleHeatmapPack=function(Data, ClusterPackResult)
{
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Re-ordered sample similarity matrix.
  
  #d=as.matrix(dist(t(Data)))
  #d=(d-min(d))/(max(d)-min(d))
  #Similarity=1-d
  Similarity=cor(Data)
  Similarity=Similarity[rownames(ClusterPackResult),rownames(ClusterPackResult)]
  return(Similarity)
}
####3.2 Pack sample similarity heatmap for already constructed simialrity matrix####
SampleHeatmapPackForSimilarity=function(Similarity, ClusterPackResult)
{
  #Input:Already constructed sample simialrity matrix.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Re-ordered sample similarity matrix.
  Similarity=Similarity[rownames(ClusterPackResult),rownames(ClusterPackResult)]
  return(Similarity)
}
####3.3 Pack sample similarity heatmap for sparse cross matrix####
SampleHeatmapPackForSparse=function(SparseCrossHeatmap)
{
  #Input:Sparse cross heatmap. Row represents clusters produced by single method, and column represents sample.
  #Output:Sample similarity matrix,which has not been re-ordered.
  Similarity=matrix(data = NA,nrow = ncol(SparseCrossHeatmap),ncol = ncol(SparseCrossHeatmap))
  colnames(Similarity)=colnames(SparseCrossHeatmap)
  rownames(Similarity)=colnames(SparseCrossHeatmap)
  for (i in 1:nrow(Similarity))
  {
    for (j in i:ncol(Similarity))
    {
      Similarity[i,j]=sum(SparseCrossHeatmap[,i]*SparseCrossHeatmap[,j],na.rm = T)
      Similarity[j,i]=Similarity[i,j]
    }
  }
  return(Similarity)
}
####4 Pack similarity network for visualization####
####4.1 Pack sample similarity network####
SampleNetworkPack=function(Similarity,ClusterPackResult)
{
  #setwd("/home/yujijun/nodejsweb/")
  #Input:Re-ordered sample similarity matrix.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Node information, edge information and filter information (for plot in webserver).
  
  #node information
  node_information=cbind(1:nrow(ClusterPackResult),rownames(ClusterPackResult),ClusterPackResult[,1])
  rownames(node_information)=NULL
  colnames(node_information)=c("id","name","cluster")
  
  #edge_information
  #source("D:/job/jobs/comsuc2.0/comsuc_pj/COMSUC数据库V1.0脚本RScript/4Methods/ConvertMatrixtoadj.R")
  diag(Similarity)=0
  Similarity[lower.tri(Similarity)]=0
  rownames(Similarity)=node_information[,1]
  colnames(Similarity)=node_information[,1]
  edge_information=ConvertMatrixtoadj(Similarity)
  edge_information[,1]=as.numeric(as.character(edge_information[,1]))
  edge_information[,2]=as.numeric(as.character(edge_information[,2]))
  index=which(edge_information[,3]==0)
  if (length(index)>0)
  {
    edge_information=edge_information[-index,]
  }
  rownames(edge_information)=NULL
  edge_information=cbind(edge_information[,1],edge_information)
  edge_information[,1]=paste(edge_information[,2],edge_information[,3],sep="_")
  colnames(edge_information)=c("id","source","target","edge_weight")
  rownames(edge_information)=NULL
  index=which(edge_information[,4]<=0)
  if (length(index)>0)
  {
    edge_information=edge_information[-index,]
  }
  rm(index)
  
  #filter information
  temp=edge_information[order(edge_information[,4],decreasing = T),]
  if (nrow(temp)>6000)
  {
    temp=temp[1:6000,]
  }
  if (min(temp[,4])<=0)
  {
    temp=temp[which(temp[,4]>0),]
  }
  filter_max=max(temp[,4])
  filter_min=min(temp[,4])
  filter_default=quantile(temp[,4],probs = 0.8)
  filter_information=matrix(data = c(filter_min,filter_default,filter_max),nrow = 1)
  colnames(filter_information)=c("edge_min","edge_default","edge_max")
  network=c(list(node_information),list(edge_information),list(filter_information))
  names(network)=c("node_information","edge_information","filter_information")
  return(network)
}
####4.2 Pack method cluster similarity network####
SampleNetworkPackForMethod=function(Similarity,Similarity_PValue,SparseCrossHeatmap,ClusterPackResult)
{
  #setwd("/home/yujijun/nodejsweb/")
  #Input:Re-ordered sample similarity matrix.
  #Input:P Value of sample similarity matrix,which has not been re-ordered.
  #Input:Sparse cross heatmap. Row represents clusters produced by single method, and column represents sample.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Node information, edge information and filter information (for plot in webserver).
  
  #node information
  node_information=cbind(1:nrow(ClusterPackResult),rownames(ClusterPackResult),ClusterPackResult[,1],rowSums(SparseCrossHeatmap[rownames(ClusterPackResult),],na.rm = T))
  rownames(node_information)=NULL
  colnames(node_information)=c("id","name","cluster","node_size")
  
  #edge_information
  #source("D:/job/jobs/comsuc2.0/comsuc_pj/COMSUC数据库V1.0脚本RScript/4Methods/ConvertMatrixtoadj.R")
  diag(Similarity)=0
  Similarity[lower.tri(Similarity)]=0
  rownames(Similarity)=node_information[,1]
  colnames(Similarity)=node_information[,1]
  edge_information=ConvertMatrixtoadj(Similarity)
  edge_information[,1]=as.numeric(as.character(edge_information[,1]))
  edge_information[,2]=as.numeric(as.character(edge_information[,2]))
  index=which(edge_information[,3]==0)
  if (length(index)>0)
  {
    edge_information=edge_information[-index,]
  }
  rownames(edge_information)=NULL
  edge_information=cbind(edge_information[,1],edge_information)
  edge_information[,1]=paste(edge_information[,2],edge_information[,3],sep="_")
  colnames(edge_information)=c("id","source","target","Jaccard_index")
  rownames(edge_information)=NULL
  edge_information=cbind(edge_information,edge_information[,4])
  colnames(edge_information)[5]="log10P"
  Similarity_PValue=Similarity_PValue[rownames(ClusterPackResult),rownames(ClusterPackResult)]
  for (i in 1:nrow(edge_information))
  {
    edge_information[i,5]=-log10(Similarity_PValue[edge_information[i,2],edge_information[i,3]])
  }
  
  #filter information
  filter_max=1
  filter_min=0
  filter_default=quantile(edge_information[,4],probs = 0.1)
  filter_information=matrix(data = c(filter_min,filter_default,filter_max),nrow = 1)
  colnames(filter_information)=c("edge_min","edge_default","edge_max")
  network=c(list(node_information),list(edge_information),list(filter_information))
  names(network)=c("node_information","edge_information","filter_information")
  return(network)
}

####5 Pack survival curve for visualization####
SurvivalPack=function(Survival_Information,ClusterPackResult)
{
  #setwd("/home/yujijun/nodejsweb/")
  #Input:Survival information matrix.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Overall pvalue, pairwise pvalue(BH ajusted pvalue) and packed survival information (for plot in webserver).
  library(survival)
  #extract useful data from survival information matrix
  Survival_Information=Survival_Information[c("_OS","_OS_IND","_EVENT","days_to_death","days_to_last_followup"),]
  Survival_Information=as.matrix(Survival_Information)
  Survival_Information["_OS_IND",which(Survival_Information["_OS_IND",]!=0)]=1
  Survival_Information["_EVENT",which(Survival_Information["_EVENT",]!=0)]=1
  n_censor=matrix(data = NA,nrow = 1,ncol = ncol(Survival_Information))
  index=which(is.na(Survival_Information["days_to_last_followup",]))
  if (length(index)>0)
  {
    n_censor[1,index]=0
  }
  index=which(!is.na(Survival_Information["days_to_last_followup",]))
  if (length(index)>0)
  {
    n_censor[1,index]=1
  }
  #pack survival information (for plot in webserver)
  Survival_Information_json=matrix(data = 0,nrow = ncol(Survival_Information),ncol = 5)
  rownames(Survival_Information_json)=colnames(Survival_Information)
  colnames(Survival_Information_json)=c("n_censor","n_event","surv","time","clust")
  Survival_Information_json[,1]=n_censor
  Survival_Information_json[,2]=t(as.numeric(Survival_Information[3,]))
  Survival_Information_json[,4]=t(as.numeric(Survival_Information[1,]))
  Survival_Information_json=Survival_Information_json[rownames(ClusterPackResult),]
  Survival_Information_json[,5]=ClusterPackResult
  index=which(is.na(Survival_Information_json[,2]))
  if (length(index)>0)
  {
    Survival_Information_json=Survival_Information_json[-index,]
  }
  index=which(Survival_Information_json[,4]==0)
  if (length(index)>0)
  {
    Survival_Information_json=Survival_Information_json[-index,]
  }
  index=which(is.na(Survival_Information_json[,4]))
  if (length(index)>0)
  {
    Survival_Information_json=Survival_Information_json[-index,]
  }
  Survival_Information_json[,"time"]=Survival_Information_json[,"time"]/30
  #compute survival ratio
  cluster_number=length(unique(Survival_Information_json[,5]))
  Survival_Information_json_list=NULL
  for (i in unique(Survival_Information_json[,5]))
  {
    index=which(Survival_Information_json[,5]==i)
    if (length(index)==1)
    {
      temp=t(as.matrix(Survival_Information_json[index,]))
      rownames(temp)=rownames(Survival_Information_json)[index]
      time_list=unique(temp[,"time"])
    }else
    {
      temp=Survival_Information_json[index,]
      temp=temp[order(temp[,4],decreasing=F),]
      time_list=unique(temp[,"time"])
    }
    S_ratio_tn=1
    for (j in 1:length(time_list))
    {
      temp1=temp[which(temp[,"time"]==time_list[j]),"n_event"]
      dn=sum(temp1,na.rm = T)#how many samples dead in the time n
      rn=length(which(temp[,"time"]>=time_list[j]))#how many samples still alive before the time n
      #rn=rn-length(which(temp[which(temp[,"time"]==time_list[j]),"n_event"]==0))#remove censor sample in time n
      S_ratio_tn1=S_ratio_tn*(1-dn/rn)
      temp[which(temp[,"time"]==time_list[j]),"surv"]=S_ratio_tn1
      S_ratio_tn=S_ratio_tn1
    }
    Survival_Information_json_list=c(Survival_Information_json_list,list(temp))
    rm(temp)
  }
  Survival_Information_json=NULL
  for (i in 1:cluster_number)
  {
    Survival_Information_json=rbind(Survival_Information_json,Survival_Information_json_list[[i]])
  }
  #compute pairwise pvalue
  library(survminer)
  res=pairwise_survdiff(Surv(time, n_event) ~ clust,data = as.data.frame(Survival_Information_json), p.adjust.method = "BH")
  pvalue_matrix=res$p.value
  #source("D:/job/jobs/comsuc2.0/comsuc_pj/COMSUC数据库V1.0脚本RScript/4Methods/ConvertMatrixtoadj.R")
  pvalue_list=ConvertMatrixtoadj(pvalue_matrix)
  index=which(is.na(pvalue_list[,3]))
  if (length(index)>0)
  {
    pvalue_list=pvalue_list[-index,]
  }
  pvalue_list=matrix(as.numeric(as.matrix(pvalue_list)),nrow = nrow(pvalue_list))
  colnames(pvalue_list)=c("i","j","p.value")
  rownames(pvalue_list)=NULL
  #compute overall pvalue
  y=Surv(time = Survival_Information_json[,4], event = Survival_Information_json[,2])
  sdf=survdiff(y ~ as.character(Survival_Information_json[,5]))
  overall_pvalue=1 - pchisq(sdf$chisq, length(sdf$n) - 1)
  
  return(list(overall_pvalue=overall_pvalue,pairwise_pvalue=pvalue_list,Survival_Information_json_list=Survival_Information_json_list))
}


####6 select signature and pack for visualization####
SignaturePack=function(Data, ClusterPackResult)
{
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Input:Cluster output packed by function ClusterPack.
  #Output:Re-ordered signature matrix.
  library(pamr)
  library(proxy)
  set.seed(120)
  Data <- Data[,rownames(ClusterPackResult)]
  mydata <- list(x=Data,y=factor(ClusterPackResult[,1]), geneid=as.character(1:nrow(Data)),
                 genenames=rownames(Data))   #reconstruction the data
  #train classifier
  mytrain<- pamr.train(mydata)
  # new.scales <- pamr.adaptthresh(mytrain)
  # mytrain2 <- pamr.train(mydata, threshold.scale=new.scales)
  
  if(!(0 %in% mytrain$nonzero))
  {
    ErrorMin = min(mytrain$errors)
    threshold <- mytrain$threshold[max(which(mytrain$errors == ErrorMin))]
    listgene <- pamr.listgenes(mytrain, mydata, threshold=threshold, genenames = T)
    signatures <- Data[listgene[,2], ]
    if(sum(mytrain$nonzero < 100) >= 1){
      UniqueError <- sort(unique(mytrain$errors))
      for(i in 1:length(UniqueError)){
        Error = UniqueError[i]
        threshold <- mytrain$threshold[max(which(mytrain$errors == Error))]
        listgene <- pamr.listgenes(mytrain, mydata, threshold=threshold, genenames = T)
        signatures <- x[listgene[,2], ]
        if(nrow(signatures) < 100){
          break
        }else{
          next
        }
      }
    }else{
      threshold <- mytrain$threshold[which.min(mytrain$nonzero)]
      listgene <- pamr.listgenes(mytrain, mydata, threshold=threshold, genenames = T)
      signatures <- Data[listgene[,2], ]
      return(signatures)
    }
  }else
  {
    ErrorMin = min(mytrain$errors)
    if(sum(mytrain$nonzero < 100) > 1){
      Less100No0_position <- setdiff(which(mytrain$nonzero < 100),which(mytrain$nonzero == 0))
      Final_position <- intersect(Less100No0_position, which(mytrain$errors == min(mytrain$errors[Less100No0_position]))) #less 100 and less error
      threshold <- mytrain$threshold[max(Final_position)]
      listgene <- pamr.listgenes(mytrain, mydata, threshold=threshold, genenames = T)
      signatures <- Data[listgene[,2], ]
    }else{
      No0_position <- mytrain$nonzero[-which(mytrain$nonzero==0)]
      threshold <- mytrain$threshold[which.min(No0_position)]
      listgene <- pamr.listgenes(mytrain, mydata, threshold=threshold, genenames = T)
      signatures <- Data[listgene[,2], ]
    }
  }
  #make sure the output of signatures without the condition of 0
  
  #sort####
  if(length(listgene[,2]) == 1)
  {
    results <- t(as.matrix(signatures))
    rownames(results)=listgene[,2]
    return(results)
  }else
  {
    #sort the samples and signatures
    d_feature = as.dist(1-cor(t(signatures),method="pearson"))
    hclust_result <- hclust(d_feature)
    signatures <- signatures[hclust_result$order, ]
    return(signatures)
  }
}

####7 Pie plot of samples for MCL####
SamplePiePack=function(ClusterPackResult)
{
  x=as.matrix(table(ClusterPackResult,useNA="ifany"))
  x=cbind(rownames(x),x)
  for (i in 1:nrow(x))
  {
    if (is.na(x[i,1]))
    {
      x[i,1]="Non-consensus"
    }else{
      x[i,1]=paste("CMS",x[i,1],sep = "")
    }
  }
  rownames(x)=NULL
  colnames(x)=c("CMS_Name","CMS_Num")
  return(x)
}

####8 Write graphML for network####
####8.1 Write graphML for method network####
Igraph_row <- function(node_row, edge_row){
  #Input: node_row and edge_row: The format is shown in the Igraph_row.RData
  #Output: A net file with igraph class
  library(igraph)
  links <- edge_row[,-1]
  nodes <- node_row
  
  #creat a net
  net <- graph_from_data_frame(
    d = links,
    vertices = nodes,
    directed = F)
  
  #set the attribute of node
  V(net)$size <- nodes$node_size
  V(net)$color <-  palette()[c(2:8,1)][nodes$cluster]
  
  #set the attribute of edge
  colors <- colorRampPalette(c("white","grey"))(101)
  E(net)$color <- colors[(floor(links$log10P)+1)]
  E(net)$width <- (links$Jaccard_index)*10
  
  #show the plot 
  # plot(
  #   net,
  #   layout=layout.fruchterman.reingold,
  #   #vertex.color = c("red","green","green","red","green","green","red"),
  #   edge.color = E(net)$color
  # )
  return(net)
}

####8.2 Write graphML for sample network####
Igraph_col <- function(node_col, edge_col, col_filter){
  #Input: node_col,edge_row and col_filter: The format is shown in the Igraph_col.RData
  #Output: A net file with igraph class
  library(igraph)
  links <- edge_col[,-1]
  links <- links[links$edge_weight >= col_filter$edge_min, ]
  nodes <- node_col
  
  #creat a net
  net <- graph_from_data_frame(
    d = links,
    vertices = nodes,
    directed = F)
  
  #set the attribute of node
  V(net)$size <- 2
  if(sum(is.na(nodes$cluster)) == 0){
    V(net)$color <-  palette()[c(2:7,1)][nodes$cluster]
    V(net)$label <- NA
  }else{
    V(net)$color <-  palette()[c(2:7,1)][nodes$cluster]
    V(net)$color[which(is.na(V(net)$color))] <- "grey"
    V(net)$label <- NA
  }
  
  #set the attribute of edge
  E(net)$width <- (links$edge_weight)*5
  
  #show the plot 
  # plot(
  #   net,
  #   #vertex.color = c("red","green","green","red","green","green","red"),
  #   layout=layout.fruchterman.reingold,
  #   node.color = V(net)$color
  # )
  return(net)
}

