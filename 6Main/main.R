##################################################################################################################################
#本地test参数
getwd()
setwd("path")
projects="TCGA"
#聚类算法选项
method=c(1,1,0,0,0,0,0,1)
names(method)=c("Kmeans", "Hclust", "SOM", "SpectralClust", "Agglomerative", "NMF", "ICA","Flame")
#指定数据
Flame_file="./Flame_BRCA_mRNA_1500.csv"
ClusterNum_input=3
cluster_file="./BRCA_mRNA_1500.txt"
clin_Information = "./clin_file.txt"
#融合算法选项
cons_method="COCA"
#读取数据
flame_cluster_out =read.table(file=Flame_file, header=TRUE, sep=",",quote="", row.names = 1,check.names = F)
FLAMENum_input =length(unique(flame_cluster_out$clusters))
Data_List=NULL
Data_List[[1]]=as.matrix(read.table(file = cluster_file ,header=T,sep="\t",quote="",check.names = F))
Data_List_back=NULL
Data_List_back[[1]]=as.matrix(read.table(file = cluster_file ,header=T,sep="\t",quote="",check.names = F))
Data_List_back[[6]]=as.matrix(read.table(file = clin_Information ,header=T,sep="\t",quote="",check.names = F))

# ##########enrich
source("./4Methods/fun_enrichgo.R")
source("./4Methods/fun_enrichkegg.R")
source("./4Methods/go_plot.R")
source("./4Methods/kegg_plot.R")

###写绝对路径
USER_DATA_GO = path + './4Methods/hsa_gene_goterm.txt'
USER_DATA_KEGG = path + './4Methods/hsa_gene_pathway.txt'
#
# ##########
start_time=Sys.time()
source("./2Platform/2.1Batch_Correction.R")
source("./3FeatureSelection/3.0FeratureSelectionFunction.R")
source("./4Methods/4.0FunctionPack.R")
source("./4Methods/methods.R")
source("./4Methods/ConvertMatrixtoadj.R")
source("./5CMS/mcl.R")
source("./5CMS/5.1SingleOmicsMultiMethod_Fusion.R")
source("./5CMS/5.1.1SingleOmicsMultiMethod_Fusion_COCA.R")
source("./5CMS/5.1.2SingleOmicsMultiMethod_Fusion_Supercluster.R")
source("./4Methods/select_sample.R")
library(igraph)

#####################################################################################################################################



if(method["Flame"]!=1)
{
  ClusterNum = ClusterNum_input
}else{
  FLAMENum = FLAMENum_input
}

#数据库和组学识别
############################
if (projects=="TCGA")
{
  omics=c(1)
  names(omics)=c("TCGA")
}
if (projects=="ICGC")
{
  omics=c(1)
  names(omics)=c("ICGC")
}
if (projects=="TARGET")
{
  omics=c(1)
  names(omics)=c("TARGET")
}
if (projects=="CPTAC")
{
  omics=c(1)
  names(omics)=c("CPTAC")
}

omics_index=which(omics==1)
method_index=which(method==1)


####标准输入数据的转置####

Data_List[[1]]=t(Data_List[[1]])

print(dim(Data_List[[1]]))



if (sum(omics)==1)
{
  if (sum(method)==1)
  {
    state_log=2
  }else
  {
    state_log=sum(method)+2
  }
}else
{
  if (sum(method)==1)
  {
    state_log=sum(omics)+2
  }else
  {
    state_log=sum(omics)*(sum(method)+1)+2
  }
}


#dir.create(Task_id)
ROOT_PATH=getwd()
write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
start_time=Sys.time()
kmeans_start_time = Sys.time()

print("Space script running name:Start Computation")
state_log_temp_progress = "Space script running progress:10|Start Computation"
print(state_log_temp_progress)
####4 Computation####
if (sum(omics)!=0)
{
  if (sum(method)==1)
  {
    ####4.1 Single-omics and Single-method####
    if (method_index==1)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=Kmeans_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)
    }
    if (method_index==2)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=Hclust_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==3)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=SOM_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==4)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=SpectralClust_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==5)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=Agglomerative_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==6)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=NMF_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==7)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=ICA_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
    }
    if (method_index==8)
    {
      state_log_temp=paste("Running ",names(method_index),sep = "")
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      result=Flame_method(Data_List[[omics_index]],cluster_result=flame_cluster_out,ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #Flame结果的加入
    }
    #output
    #setwd(paste("./", Task_id,sep = ""))
    dir.create(names(omics_index))
    setwd(paste("./", names(omics_index),sep = ""))
    dir.create(names(method_index))
    setwd(paste("./", names(method_index),sep = ""))
    dir.create("heatmap")
    dir.create("network")
    dir.create("survival")
    dir.create("biomarker")
    dir.create("enrichment")
    #4.1.1 Clusters and Heatmap
    temp=result$cluster_output
    temp=cbind(rownames(temp),temp)
    colnames(temp)=c("name","cluster")
    write.table(temp,file = paste("heatmap/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
    temp2=result$heatmap_result
    b = split(1:dim(temp)[1], temp[,2])
    for (i in 1:length(b)){
      if (length(b[[i]]) == 1)
      {
        print("cluster number is 1")
      }
      else
      {
        c_node = as.matrix(temp[b[[i]],2])
        #富集分析
        genes = row.names(c_node)
        #聚类数id
        clusterid = c_node[1]
        c=cbind(rownames(c_node),c_node)
        colnames(c)=c("name","cluster")
        dir.create(paste("heatmap","/heatmapCluster",clusterid,sep=""))
        write.table(c,file = paste("heatmap/heatmapCluster",clusterid,"/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
        d = as.array(rownames(c))
        e = temp2[d,d]
        write.table(e,file = paste("heatmap/heatmapCluster",clusterid,"/heatmap.txt",sep=""),quote = F,sep = "\t")
        network_result_cluster = SampleNetworkPack(Similarity=e,ClusterPackResult=c_node)
        ########Survival
        survival_result = Survival(Data_List_back[[1]],Data_List_back[[6]],network_result_cluster$node_information)
        survival_result_modue = Survival_modue(Data_List_back[[1]],Data_List_back[[6]],network_result_cluster$node_information)
        #获得基因名
        Survival_result_gene_name = Survival_gene_name(Data_List_back[[1]],network_result_cluster$node_information)
        #处理特殊字符
        Survival_result_gene_name = lapply(Survival_result_gene_name, function(x) gsub("[&!#?|]", "", x))
        dir.create(paste("survival","/survivalCluster",clusterid,sep=""))
        genes_pvalue_ranks = NULL
        for (g in 1:length(Survival_result_gene_name)){
          dir.create(paste("survival","/survivalCluster", clusterid, "/", Survival_result_gene_name[[g]],sep=""))
          survival_result_gene = survival_result[[g]]
          survival_result_modue_gene = survival_result_modue[[g]]
          genes_pvalue_ranks[g] = paste(Survival_result_gene_name[g],survival_result_gene$overall_pvalue,survival_result_modue_gene$overall_pvalue, sep=",")
          write.table(survival_result_gene$overall_pvalue,file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
          write.table(survival_result_modue_gene$overall_pvalue,file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
          write.table(survival_result_gene$pairwise_pvalue,file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
          write.table(survival_result_modue_gene$pairwise_pvalue,file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
          for (k in 1:length(survival_result_gene$Survival_Information_json_list))
          {
            write.table(survival_result_gene$Survival_Information_json_list[[k]],file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/survival_",k,".txt",sep=""),quote = F,sep = "\t",row.names = F)
            write.table(survival_result_modue_gene$Survival_Information_json_list[[k]],file = paste("survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_survival_",k,".txt",sep=""),quote = F,sep = "\t",row.names = F)
          }
        }
        write.table(genes_pvalue_ranks,file = paste("survival/survivalCluster",clusterid,"/genes_pvalue_ranks2.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
        ########
        dir.create(paste("network","/networkCluster",clusterid,sep=""))
        write.table(network_result_cluster$node_information,file = paste("network/networkCluster",clusterid,"/node.txt",sep=""),quote = F,sep = "\t",row.names = F)
        write.table(network_result_cluster$edge_information,file = paste("network/networkCluster",clusterid,"/edge.txt",sep=""),quote = F,sep = "\t",row.names = F)
        write.table(network_result_cluster$filter_information,file = paste("network/networkCluster",clusterid,"/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
        net_cluster=Igraph_col(node_col=as.data.frame(network_result_cluster$node_information),edge_col=network_result_cluster$edge_information, col_filter=as.data.frame(network_result_cluster$filter_information))
        write_graph(net_cluster, file = paste("network/networkCluster",clusterid,"/network.graphml",sep=""), format = c("graphml"))
        #################################enrich
        library(AnnotationDbi)
        library(org.Hs.eg.db)
        mapid = mapIds(org.Hs.eg.db,genes,"ENTREZID","SYMBOL")
        mapid = na.omit(mapid)
        source("D:/job/jobs/comsuc2.0/run_comsuc2/4Methods/fun_enrichkegg.R")
        source("D:/job/jobs/comsuc2.0/run_comsuc2/4Methods/fun_enrichgo.R")
        if(length(mapid)>0){
          kegg = fun_enrich_internal(gene = mapid,USER_DATA = USER_DATA_KEGG)
          # KEGG
          if(nrow(kegg) > 0){
            dir.create(paste("enrichment","/KeggCluster",clusterid,sep=""))
            write.table(kegg,file = paste("enrichment/KeggCluster",clusterid,"/Table_KEGG.csv",sep=""),quote = F,sep = ",",col.names = T,row.names = F)
            keggplot = plot_keggdot(data = kegg)
            ggsave(keggplot,filename = paste("enrichment/KeggCluster",clusterid,"/kegg_barplot.png",sep=""),dpi = 600, height = 8)
          }else{
            message("--> No significant KEGG pathway....")
          }
          
          # GO
          go = fun_enrich_go(gene = mapid,ont = 'ALL',USER_DATA = USER_DATA_GO)
          # GO输出结果
          if(nrow(go) > 0){
            dir.create(paste("enrichment","/GoCluster",clusterid,sep=""))
            write.table(go,file = paste("enrichment/GoCluster",clusterid,"/Table_GO.csv",sep=""),quote = F,sep = ",",col.names = T,row.names = F)
            goplot = gobarplot(data = go)
            ggsave(goplot,filename = paste("enrichment/GoCluster",clusterid,"/go_barplot.png",sep=""),dpi = 600, width = 8, height = 9)
          }else{
            message("--> fun_enrich_go: No significat GO term....")
          }
        }else{
          message("--> fun_enrich_go: Missing gene parameter....")
        }
        rm(genes)
        rm(mapid)
        rm(c)
        rm(e)
        rm(c_node)
        rm(network_result_cluster)
        rm(net_cluster)
      }
    }
    rm(b)
    rm(temp)
    write.table(result$heatmap_result,file = paste("heatmap/heatmap.txt",sep=""),quote = F,sep = "\t")
    # #4.1.2 Network
    write.table(result$network_result$node_information,file = paste("network/node.txt",sep=""),quote = F,sep = "\t",row.names = F)
    write.table(result$network_result$edge_information,file = paste("network/edge.txt",sep=""),quote = F,sep = "\t",row.names = F)
    write.table(result$network_result$filter_information,file = paste("network/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
    net=Igraph_col(node_col=as.data.frame(result$network_result$node_information),edge_col=result$network_result$edge_information, col_filter=as.data.frame(result$network_result$filter_information))
    write_graph(net, file = paste("network/network.graphml",sep=""), format = c("graphml"))
    rm(net)
    #4.1.4 Signature
    write.table(result$signature_result,file = paste("biomarker/biomarker.txt",sep=""),quote = F,sep = "\t")
    temp=result$cluster_output
    temp=cbind(rownames(temp),temp)
    colnames(temp)=c("name","cluster")
    write.table(temp,file = paste("biomarker/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
    rm(temp)
    end_time=Sys.time()
    all_time=end_time-start_time
    #save(Data_List,result,cancer_type,method,omics,projects,start_time,end_time,all_time,file = paste(ROOT_PATH,"/Result.RData",sep = ""))
    state_log=rbind(state_log,"End")
    write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  }else
  {
    ####4.2 Single-omics and Multi-method integration####
    res_single_list=NULL
    if (length(which(method_index==1)))
    {
      state_log_temp="Running Kmeans"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=Kmeans_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(Kmeans=temp))
      rm(temp)
    }
    if (length(which(method_index==2)))
    {
      state_log_temp="Running iKmeans"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=Hclust_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(Hclust=temp))
      rm(temp)
    }
    if (length(which(method_index==3)))
    {
      state_log_temp="Running Hclust"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=SOM_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(SOM=temp))
      rm(temp)
    }
    if (length(which(method_index==4)))
    {
      state_log_temp="Running iHclust"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=SpectralClust_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(SpectralClust=temp))
      rm(temp)
    }
    if (length(which(method_index==5)))
    {
      state_log_temp="Running NMF"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=Agglomerative_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(Agglomerative=temp))
      rm(temp)
    }
    if (length(which(method_index==6)))
    {
      state_log_temp="Running SOM"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=NMF_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(NMF=temp))
      rm(temp)
    }
    if (length(which(method_index==7)))
    {
      state_log_temp="Running Spectralclust"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=ICA_method_v2(Data_List[[omics_index]], ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #输入转置后的矩阵，对gene聚类
      res_single_list=c(res_single_list,list(ICA=temp))
      rm(temp)
    }
    if (length(which(method_index==8)))
    {
      state_log_temp="Running Spectralclust"
      state_log=rbind(state_log,state_log_temp)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(state_log_temp)
      temp=Flame_method(Data_List[[omics_index]],cluster_result=flame_cluster_out,ClusterNum = ClusterNum, FLAMENum = FLAMENum)  #Flame结果的加入
      res_single_list=c(res_single_list,list(Flame=temp))
      rm(temp)
    }
    kmeans_end_time = Sys.time()
    kmeans_all_time = kmeans_end_time-kmeans_start_time
    print("算法时间")
    print(kmeans_all_time)
    state_log_temp="Running Multi-methods Integration"
    state_log_time="Multi-methods Analysis time"
    state_log_temp_progress = "Space script running progress:30|Consensus Analysis"
    print(state_log_temp_progress)
    state_log=rbind(state_log,state_log_temp,state_log_time,kmeans_all_time,state_log_temp_progress)
    write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
    rm(state_log_temp)
    #######MCL
    mcl_start_time = Sys.time()
    if (cons_method=="MCL")
    {
      MCL_result=SOMMFusionPackv2(ResultSingleMethod=res_single_list, Data=Data_List[[omics_index]], ClustMethod=names(res_single_list), OmicsLevel=names(omics_index))
    }
    if (cons_method=="COCA")
    {
      MCL_result=SOMMFusionPack_COCAv2(ResultSingleMethod=res_single_list, Data=Data_List[[omics_index]], ClustMethod=names(res_single_list), OmicsLevel=names(omics_index))
    }
    if (cons_method=="Supercluster")
    {
      MCL_result=SOMMFusionPack_Superclusterv2(ResultSingleMethod=res_single_list, Data=Data_List[[omics_index]], ClustMethod=names(res_single_list), OmicsLevel=names(omics_index))
    }
    mcl_end_time = Sys.time()
    mcl_time = mcl_end_time - mcl_start_time
    print("MCL算法时间")
    print(mcl_time)
    state_log_time="Cons Analysis"
    state_log_temp_progress = "Space script running progress:60|Gene and Survival Analysis"
    print(state_log_temp_progress)
    state_log=rbind(state_log,state_log_time,mcl_time,state_log_temp_progress)
    write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
    #output
    #setwd(paste("./", Task_id,sep = ""))
    ####4.2.1 Single method result output#####
    dir.create(names(omics_index))
    setwd(paste("./", names(omics_index),sep = ""))
    m=1
    for (j in 1:length(method))
    {
      if (method[j]!=0)
      {
        dir.create(names(method)[j])
        dir.create(paste(names(method)[j],"/heatmap",sep=""))
        dir.create(paste(names(method)[j],"/network",sep=""))
        dir.create(paste(names(method)[j],"/survival",sep=""))
        dir.create(paste(names(method)[j],"/biomarker",sep=""))
        dir.create(paste(names(method)[j],"/enrichment",sep=""))
        result=res_single_list[[m]]
        #Clusters and Heatmap
        temp=result$cluster_output
        temp=cbind(rownames(temp),temp)
        colnames(temp)=c("name","cluster")
        write.table(temp,file = paste(names(method)[j],"/heatmap/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
        temp2=result$heatmap_result
        rm(temp)
        rm(temp2)
        write.table(result$heatmap_result,file = paste(names(method)[j],"/heatmap/heatmap.txt",sep=""),quote = F,sep = "\t")
        #Network
        write.table(result$network_result$node_information,file = paste(names(method)[j],"/network/node.txt",sep=""),quote = F,sep = "\t",row.names = F)
        write.table(result$network_result$edge_information,file = paste(names(method)[j],"/network/edge.txt",sep=""),quote = F,sep = "\t",row.names = F)
        write.table(result$network_result$filter_information,file = paste(names(method)[j],"/network/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
        net=Igraph_col(node_col=as.data.frame(result$network_result$node_information),edge_col=result$network_result$edge_information, col_filter=as.data.frame(result$network_result$filter_information))
        write_graph(net, file = paste(names(method)[j],"/network/network.graphml",sep=""), format = c("graphml"))
        rm(net)
        #Signature
        write.table(result$signature_result,file = paste(names(method)[j],"/biomarker/biomarker.txt",sep=""),quote = F,sep = "\t")
        temp=result$cluster_output
        temp=cbind(rownames(temp),temp)
        colnames(temp)=c("name","cluster")
        write.table(temp,file = paste(names(method)[j],"/biomarker/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
        rm(temp)
        rm(result)
        m=m+1
      }
    }
    rm(m)
    ####4.2.2 consensus result output based on multi-method####
    dir.create("consensus")
    if (is.atomic(MCL_result))
    {
      write.table(MCL_result,file = paste("consensus/error.txt",sep = ""),quote = F,col.names = F,row.names = F,sep = "\t")
    }else
    {
      dir.create("consensus/row_heatmap")
      dir.create("consensus/row_network")
      dir.create("consensus/col_heatmap")
      dir.create("consensus/col_network")
      dir.create("consensus/cross_heatmap")
      dir.create("consensus/pie")
      dir.create("consensus/survival")
      dir.create("consensus/biomarker")
      dir.create("consensus/enrichment")
      ####4.2.2.1 Method result####
      temp=MCL_result$method_list$method_cluster_output
      temp=cbind(rownames(temp),temp)
      colnames(temp)=c("name","cluster")
      write.table(temp,file = paste("consensus/row_heatmap/clusters.txt",sep=""),quote = F,sep = "\t",row.names = F)
      rm(temp)
      write.table(MCL_result$method_list$method_heatmap_result,file = paste("consensus/row_heatmap/heatmap_row.txt",sep=""),quote = F,sep = "\t",row.names = T)
      write.table(MCL_result$method_list$method_information,file = paste("consensus/row_heatmap/NodeInformation.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$method_list$method_network_result$node_information,file = paste("consensus/row_network/node_row.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$method_list$method_network_result$edge_information,file = paste("consensus/row_network/edge_row.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$method_list$method_network_result$filter_information,file = paste("consensus/row_network/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
      net=Igraph_row(node_row=as.data.frame(MCL_result$method_list$method_network_result$node_information),edge_row=MCL_result$method_list$method_network_result$edge_information)
      write_graph(net, file = paste("consensus/row_network/network.graphml",sep=""), format = c("graphml"))
      rm(net)
      ####4.2.2.2 Sample result####
      temp=MCL_result$sample_list$sample_cluster_output
      temp=cbind(rownames(temp),temp)
      colnames(temp)=c("name","cluster")
      write.table(temp,file = paste("consensus/col_heatmap/clusters.txt",sep=""),quote = F,sep = "\t",row.names = F)
      cluster_num = split(1:dim(temp)[1], temp[,2])
      temp2 = MCL_result$sample_list$sample_heatmap_result
      analysis_start_time = Sys.time()
      for (i in 1:length(cluster_num)){
        if (length(cluster_num[[i]]) == 1)
        {
          print("cluster number is 1")
        }
        else
        {
          c_node = as.matrix(temp[cluster_num[[i]],2])
          genes = row.names(c_node)
          clusterid = c_node[1]
          c=cbind(rownames(c_node),c_node)
          colnames(c)=c("name","cluster")
          dir.create(paste("consensus/col_heatmap","/heatmapCluster",clusterid,sep=""))
          write.table(c,file = paste("consensus/col_heatmap","/heatmapCluster",clusterid,"/clusters.txt",sep=""),quote = F,sep = "\t",col.names = T,row.names = F)
          d = as.array(rownames(c))
          e = temp2[d,d]
          write.table(e,file = paste("consensus/col_heatmap","/heatmapCluster",clusterid,"/heatmap_col.txt",sep=""),quote = F,sep = "\t")
          network_cluster = SampleNetworkPack(Similarity=e,ClusterPackResult=c_node)
          ########Survival
          surv_start_time = Sys.time()
          survival_result = Survival(Data_List_back[[1]],Data_List_back[[6]],network_cluster$node_information)
          survival_result_modue = Survival_modue(Data_List_back[[1]],Data_List_back[[6]],network_cluster$node_information)
          #获得基因名
          Survival_result_gene_name = Survival_gene_name(Data_List_back[[1]],network_cluster$node_information)
          #处理特殊字符
          Survival_result_gene_name = lapply(Survival_result_gene_name, function(x) gsub("[&!#?|]", "", x))
          dir.create(paste("consensus/survival","/survivalCluster",clusterid,sep=""))
          #获得基因的pvalue_rank
          genes_pvalue_ranks = NULL
          for (g in 1:length(Survival_result_gene_name)){
            dir.create(paste("consensus/survival","/survivalCluster", clusterid, "/", Survival_result_gene_name[[g]],sep=""))
            survival_result_gene = survival_result[[g]]
            survival_result_modue_gene = survival_result_modue[[g]]
            #获得基因的pvalue_rank
            genes_pvalue_ranks[g] = paste(Survival_result_gene_name[g],survival_result_gene$overall_pvalue,survival_result_modue_gene$overall_pvalue, sep=",")
            #
            write.table(survival_result_gene$overall_pvalue,file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
            write.table(survival_result_modue_gene$overall_pvalue,file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
            write.table(survival_result_gene$pairwise_pvalue,file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
            write.table(survival_result_modue_gene$pairwise_pvalue,file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
            for (k in 1:length(survival_result_gene$Survival_Information_json_list))
            {
              write.table(survival_result_gene$Survival_Information_json_list[[k]],file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/survival_",k,".txt",sep=""),quote = F,sep = "\t",row.names = F)
              write.table(survival_result_modue_gene$Survival_Information_json_list[[k]],file = paste("consensus/survival/survivalCluster",clusterid,"/",Survival_result_gene_name[[g]],"/modue_survival_",k,".txt",sep=""),quote = F,sep = "\t",row.names = F)
            }
          }
          #保存的pvalue_rank
          write.table(genes_pvalue_ranks,file = paste("consensus/survival/survivalCluster",clusterid,"/genes_pvalue_ranks.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
          ########
          #surv_start_time = Sys.time()
          surv_end_time = Sys.time()
          surv_time = surv_end_time - surv_start_time
          print("surv算法时间")
          print(surv_time)
          state_log_temp="Running Survival Analysis"
          state_log_time="Survival Analysis time"
          state_log=rbind(state_log,state_log_temp,state_log_time,surv_time)
          write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
          rm(state_log_temp)
          ########
          ########
          dir.create(paste("consensus/col_network","/networkCluster",clusterid,sep=""))
          write.table(network_cluster$node_information,file = paste("consensus/col_network/networkCluster",clusterid,"/node_col.txt",sep=""),quote = F,sep = "\t",row.names = F)
          write.table(network_cluster$edge_information,file = paste("consensus/col_network/networkCluster",clusterid,"/edge_col.txt",sep=""),quote = F,sep = "\t",row.names = F)
          write.table(network_cluster$filter_information,file = paste("consensus/col_network/networkCluster",clusterid,"/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
          net_cluster=Igraph_col(node_col=as.data.frame(network_cluster$node_information),edge_col=network_cluster$edge_information, col_filter=as.data.frame(network_cluster$filter_information))
          write_graph(net_cluster, file = paste("consensus/col_network/networkCluster",clusterid,"/network.graphml",sep=""), format = c("graphml"))
          #################################enrich
          enrich_start_time = Sys.time()
          library(AnnotationDbi)
          library(org.Hs.eg.db)
          mapid = mapIds(org.Hs.eg.db,genes,"ENTREZID","SYMBOL")
          mapid = na.omit(mapid)
          if(length(mapid)>0){
            kegg = fun_enrich_internal(gene = mapid,USER_DATA = USER_DATA_KEGG)
            # KEGG
            if(nrow(kegg) > 0){
              dir.create(paste("consensus/enrichment","/KeggCluster",clusterid,sep=""))
              #write.csv(kegg,"Table_KEGG.csv",row.names = F)
              write.table(kegg,file = paste("consensus/enrichment/KeggCluster",clusterid,"/Table_KEGG.csv",sep=""),quote = F,sep = ",",col.names = T,row.names = F)
              keggplot = plot_keggdot(data = kegg)
              ggsave(keggplot,filename = paste("consensus/enrichment/KeggCluster",clusterid,"/kegg_barplot.png",sep=""),dpi = 600, height = 8)
            }else{
              message("--> No significant KEGG pathway....")
            }
            # GO
            go = fun_enrich_go(gene = mapid,ont = 'ALL',USER_DATA = USER_DATA_GO)
            # GO输出结果
            if(nrow(go) > 0){
              dir.create(paste("consensus/enrichment","/GoCluster",clusterid,sep=""))
              write.table(go,file = paste("consensus/enrichment/GoCluster",clusterid,"/Table_GO.csv",sep=""),quote = F,sep = ",",col.names = T,row.names = F)
              goplot = gobarplot(data = go)
              ggsave(goplot,filename = paste("consensus/enrichment/GoCluster",clusterid,"/go_barplot.png",sep=""),dpi = 600, width = 8, height = 9)
            }else{
              message("--> fun_enrich_go: No significat GO term....")
            }
          }
          #enrich_start_time = Sys.time()
          enrich_end_time = Sys.time()
          enrich_time = enrich_end_time - enrich_start_time
          print("enrich算法时间")
          print(enrich_time)
          state_log_temp="Running Enrich Analysis"
          state_log_time="Enrich Analysis time"
          state_log=rbind(state_log,state_log_temp,state_log_time,enrich_time)
          write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
          rm(state_log_temp)
          rm(genes)
          rm(mapid)
          rm(c)
          rm(e)
          rm(c_node)
          rm(network_cluster)
          rm(net_cluster)
        }
      }
      #analysis_start_time = Sys.time()
      analysis_end_time = Sys.time()
      analysis_time = analysis_end_time - analysis_start_time
      print("analysis算法时间")
      print(analysis_time)
      state_log_temp="Running Analysis"
      state_log_time="Analysis time"
      state_log_temp_progress = "Space script running progress:90|Result Saving"
      print(state_log_temp_progress)
      state_log=rbind(state_log,state_log_temp,state_log_time,analysis_time,state_log_temp_progress)
      write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
      rm(temp)
      rm(temp2)
      write.table(MCL_result$sample_list$sample_heatmap_result,file = paste("consensus/col_heatmap/heatmap_col.txt",sep=""),quote = F,sep = "\t",row.names = T)
      write.table(MCL_result$sample_list$sample_network_result$node_information,file = paste("consensus/col_network/node_col.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$sample_list$sample_network_result$edge_information,file = paste("consensus/col_network/edge_col.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$sample_list$sample_network_result$filter_information,file = paste("consensus/col_network/filter.txt",sep=""),quote = F,sep = "\t",row.names = F)
      write.table(MCL_result$sample_list$sample_pie_result,file = paste("consensus/pie/pie.txt",sep=""),quote = F,sep = "\t",row.names = F)
      net=Igraph_col(node_col=as.data.frame(MCL_result$sample_list$sample_network_result$node_information),edge_col=MCL_result$sample_list$sample_network_result$edge_information, col_filter=as.data.frame(MCL_result$sample_list$sample_network_result$filter_information))
      write_graph(net, file = paste("consensus/col_network/network.graphml",sep=""), format = c("graphml"))
      rm(net)
      ####4.2.2.3 Cross heatmap####
      write.table(MCL_result$cross_heatmap_list$cross_heatmap,file = paste("consensus/cross_heatmap/cross_heatmap.txt",sep=""),quote = F,sep = "\t",row.names = T)
      write.table(MCL_result$cross_heatmap_list$sparse_cross_heatmap,file = paste("consensus/cross_heatmap/sparse_cross_heatmap.txt",sep=""),quote = F,sep = "\t",row.names = T)
      temp=MCL_result$sample_list$sample_cluster_output
      temp=cbind(rownames(temp),temp)
      colnames(temp)=c("name","cluster")
      write.table(temp,file = paste("consensus/cross_heatmap/clusters.txt",sep=""),quote = F,sep = "\t",row.names = F)
      rm(temp)
      ####4.2.2.5 Signature####
      write.table(MCL_result$signature_result,file = paste("consensus/biomarker/biomarker.txt",sep=""),quote = F,sep = "\t")
      temp=MCL_result$sample_list$sample_cluster_output
      temp=cbind(rownames(temp),temp)
      colnames(temp)=c("name","cluster")
      index_temp=which(is.na(temp[,2]))
      if (length(index_temp)>0)
      {
        temp=temp[-index_temp,]
      }
      rm(index_temp)
      write.table(temp,file = paste("consensus/biomarker/clusters.txt",sep=""),quote = F,sep = "\t",row.names = F)
      rm(temp)
    }
    end_time=Sys.time()
    all_time=end_time-start_time
    print("总的时间")
    print(all_time)
    state_log_temp_progress = "Space script running progress:100|End"
    print(state_log_temp_progress)
    save(Data_List,MCL_result,res_single_list,method,omics,projects,start_time,end_time,all_time,file = paste(ROOT_PATH,"/Result.RData",sep = ""))
    state_log=rbind(state_log,"End",all_time,state_log_temp_progress)
    write.table(state_log,file = paste(ROOT_PATH,"/Log.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = F)
  }
}

