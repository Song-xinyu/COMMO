SOMMFusionPack_Superclusterv2=function(ResultSingleMethod, Data, ClustMethod, OmicsLevel){
  #This function is used for integrating results computed by each method (single omics data).
  #Input:Result list produced by each single method.
  #Input:Clinical Information of each sample.
  #Input:Data matrix. Row represents feature, and column represents samples.
  #Input:Each method name.
  #Input:Omics Data name.
  #Output:Pack of results.
  #setwd("/home/lihl/comsuc_api/")
  library(proxy)
  #library(MCL)
  library(WeightedCluster)
  #library(stringr)
  #library(survival)
  # source("./4Methods/4.0FunctionPack.R")
  # source("./5CMS/mcl.R")
  ####1 Cross heatmap matrix for each result clustering based on single method####
  sparse_cross_heatmap=NULL
  cross_heatmap=NULL
  for (i in 1:length(ClustMethod))
  {
    new_matrix=matrix(0, ResultSingleMethod[[i]]$cluster_number, ncol(Data))
    new_matrix_temp=matrix(0,1,ncol(Data))
    colnames(new_matrix)=colnames(Data)
    colnames(new_matrix_temp)=colnames(Data)
    rownames(new_matrix_temp)=ClustMethod[i]
    rownames(new_matrix)=paste0(OmicsLevel,"_",ClustMethod[i], "_", c(1:ResultSingleMethod[[i]]$cluster_number))
    for(j in 1:ResultSingleMethod[[i]]$cluster_number)
    {
      new_matrix[j,rownames(ResultSingleMethod[[i]]$cluster_output)[which(ResultSingleMethod[[i]]$cluster_output==j)]]=1
      new_matrix_temp[1,rownames(ResultSingleMethod[[i]]$cluster_output)[which(ResultSingleMethod[[i]]$cluster_output==j)]]=j
    }
    sparse_cross_heatmap=rbind(sparse_cross_heatmap, new_matrix)
    cross_heatmap=rbind(cross_heatmap,new_matrix_temp)
    rm(new_matrix)
    rm(new_matrix_temp)
  }
  
  ####2 Jaccard similarity matrix and adjust p-value matrix for each result clustering based on single method####
  PVALUE_CUTOFF=0.05
  cluster_jaccard=1-as.matrix(dist(sparse_cross_heatmap, method = "Jaccard"))
  cluster_jaccard_pvalue=cluster_jaccard
  cluster_jaccard_pvalue[which(cluster_jaccard_pvalue!=0)]=NA
  for (i in 1:nrow(cluster_jaccard_pvalue))
  {
    for (j in 1:ncol(cluster_jaccard_pvalue))
    {
      M=sum(sparse_cross_heatmap[i,])
      n=sum(sparse_cross_heatmap[j,])
      N=ncol(sparse_cross_heatmap)
      k=length(which(sparse_cross_heatmap[i,]*sparse_cross_heatmap[j,]==1))
      cluster_jaccard_pvalue[i,j]=phyper(k-1,M, N-M, n, lower.tail=FALSE)
      rm(M,N,n,k)
    }
  }
  cluster_jaccard_pvalue=matrix(p.adjust(cluster_jaccard_pvalue,method = "BH"),nrow = nrow(cluster_jaccard))
  rownames(cluster_jaccard_pvalue)=rownames(cluster_jaccard)
  colnames(cluster_jaccard_pvalue)=colnames(cluster_jaccard)
  cluster_jaccard[which(cluster_jaccard_pvalue>PVALUE_CUTOFF)]=0#Filter the Jaccard similarity edge whose pvalue is higher than cutoff
  diag(cluster_jaccard)=1
  #pheatmap(cluster_jaccard)
  
  #####bug#####  修复去重标准差为0的行列    
  data = sparse_cross_heatmap
  sparse_cross_heatmap <- data[apply(data, 1, function(x) sd(x)!=0),]
  ######
  ####3 Supercluster####
  if (nrow(sparse_cross_heatmap)<=20)
  {
    cluster_number=2
  }else
  {
    cluster_number=OptClustNum(Data=t(sparse_cross_heatmap),Method="hc") #Determine the optimal cluster number
  }
  hc=hclust(as.dist(1-cor(t(sparse_cross_heatmap))),"ave") #Using the function of kmeans to calculate the cluster
  cluster_result=unname(cutree(hc,k=cluster_number))
  rm(hc)
  
  ####5 method cluster output####
  cluster_output=ClusterPack(Data=t(sparse_cross_heatmap),cluster=cluster_result)
  ####6 method heatmap####
  heatmap_result=SampleHeatmapPackForSimilarity(Similarity=cluster_jaccard,ClusterPackResult=cluster_output)
  #pheatmap(heatmap_result,cluster_rows = F,cluster_cols = F,annotation_col=as.data.frame(cluster_output))
  ####7 method network####
  network_result=SampleNetworkPackForMethod(Similarity=heatmap_result,Similarity_PValue=cluster_jaccard_pvalue,SparseCrossHeatmap=sparse_cross_heatmap,ClusterPackResult=cluster_output)
  ####8 sample similarity heatmap####
  sample_heatmap_result=SampleHeatmapPackForSparse(SparseCrossHeatmap=sparse_cross_heatmap)/nrow(cross_heatmap)
  #pheatmap(sample_heatmap_result,cluster_rows = F,cluster_cols = F)
  ####9 identify sample CMS####
  if (length(ClustMethod)<=3)
  {
    sample_cluster_output=SampleCMSPackFortwothree(SparseCrossHeatmap=sparse_cross_heatmap,ClusterPackResult=cluster_output,ClustMethod=ClustMethod)
  }else
  {
    sample_cluster_output=SampleCMSPack(SparseCrossHeatmap=sparse_cross_heatmap,ClusterPackResult=cluster_output,ClustMethod=ClustMethod)
  }
  ########
  if (length(which(is.na(sample_cluster_output)))==nrow(sample_cluster_output))
  {
    return("Error! All non-consensus!")
  }else{
    if (length(which(is.na(sample_cluster_output)))>0)
    {
      num_temp=length(unique(sample_cluster_output))-1
    }else{
      num_temp=length(unique(sample_cluster_output))
    }
    if (num_temp==1)
    {
      return("Error! Only one CMS!")
    }else{
    ####10 re-ordered sample similarity heatmap####
    sample_heatmap_result=SampleHeatmapPackForSimilarity(Similarity=sample_heatmap_result, ClusterPackResult=sample_cluster_output)
    #pheatmap(sample_heatmap_result,cluster_rows = F,cluster_cols = F,annotation_col=as.data.frame(sample_cluster_output))
    ####11 sample similarity network####
    sample_network_result=SampleNetworkPack(Similarity=sample_heatmap_result,ClusterPackResult=sample_cluster_output)
    ####12 pie plot of samples####
    sample_pie_result=SamplePiePack(ClusterPackResult=sample_cluster_output)
    ####13 Pack Cross heatmap####
    sparse_cross_heatmap=sparse_cross_heatmap[rownames(cluster_output),rownames(sample_cluster_output)]
    cross_heatmap=cross_heatmap[,rownames(sample_cluster_output)]
    #pheatmap(cross_heatmap,cluster_rows = F,cluster_cols = F,annotation_col=as.data.frame(sample_cluster_output))
    #pheatmap(sparse_cross_heatmap,cluster_rows = F,cluster_cols = F,annotation_col=as.data.frame(sample_cluster_output),annotation_row=as.data.frame(cluster_output))
    ####14 Information####
    ####14.1 Method Information####
    Method_Information=matrix(data = NA,nrow = nrow(cluster_output),ncol = 5)
    colnames(Method_Information)=c("Node Name","Consensue Subtypes","Cluster Method","Omics Data","Sample Number")
    rownames(Method_Information)=NULL
    Method_Information[,c(1,2,5)]=network_result$node_information[,2:4]
    Method_Information[,4]=OmicsLevel
    for (i in 1:nrow(Method_Information))
    {
      if (is.na(Method_Information[i,2]))
      {
        Method_Information[i,2]="Non-consensus"
      }else
      {
        Method_Information[i,2]=paste("CMS",Method_Information[i,2],sep = "")
      }
      Method_Information[i,3]=unlist(strsplit(Method_Information[i,1], "_"))[2]
    }
    ####14.2 Sample Information####
    # Sample_Information=SurvivalInformation
    # Sample_Information=Sample_Information[c("_cohort","gender","age_at_initial_pathologic_diagnosis","days_to_last_followup",
    #                                         "days_to_death","vital_status","clinical_M","laterality"),rownames(sample_cluster_output)]
    # rownames(Sample_Information)=c("Cohort","Gender","Age at Initial Diagnosis","Days to Last Followup",
    #                                "Days to Death","Vital Status","Clinical M","Laterality")
    ####15 Survival####
    # na_index=which(is.na(sample_cluster_output))
    # if (length(na_index)>0)
    # {
    #   sample_cluster_output_temp=as.matrix(sample_cluster_output[-na_index,])
    # }else
    # {
    #   sample_cluster_output_temp=sample_cluster_output
    # }
    # survival_result=SurvivalPack(Survival_Information=SurvivalInformation,ClusterPackResult=sample_cluster_output_temp)
    ####16 Signature####
    #sample_cluster_output_temp=sample_cluster_output
    #if (length(na_index)>0)
    #{
    #  sample_cluster_output_temp[na_index,]="Non-Consensus"
    #}
    # signature_result=SignaturePack(Data=Data,ClusterPackResult=sample_cluster_output_temp)
    #pheatmap(signature_result,cluster_rows = F,cluster_cols = F,annotation_col=as.data.frame(sample_cluster_output_temp))
    # rm(sample_cluster_output_temp)
    
    method_list=list(method_cluster_output=cluster_output,method_heatmap_result=heatmap_result,
                     method_network_result=network_result,method_information=Method_Information)
    sample_list=list(sample_cluster_output=sample_cluster_output,sample_heatmap_result=sample_heatmap_result,
                     sample_network_result=sample_network_result,sample_pie_result=sample_pie_result)
    cross_heatmap_list=list(cross_heatmap=cross_heatmap,sparse_cross_heatmap=sparse_cross_heatmap)
    result_all=list(method_list=method_list,sample_list=sample_list,cross_heatmap_list=cross_heatmap_list)
    return(result_all)
    }
  }
}




