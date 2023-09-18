# setwd('D:/job/jobs/comsuc2.0/run_comsuc2/预后模型/select-sample')
#生存单基因表达算法
Survival = function(Data,Survival_Information,node_col){
  library(readxl)
  Data = Data
  Survival_Information = Survival_Information
  node_col = node_col
  df1 = as.data.frame(Data)
  df2 = as.data.frame(node_col)
  df <- df1[rownames(df1) %in% df2$name,]
  #1.按照单基因排序，筛选出大于上四分位数和小于下四分位数的样本
  survival_result = NULL
  for (gene in 1:nrow(df)){
    df_gene <-df[gene,]
    # 计算标准差的上四分位数和下四分位数
    upper_quartile <- apply(df_gene, 1, quantile, probs = 0.75)
    lower_quartile <- apply(df_gene, 1, quantile, probs = 0.25)
    # 筛选出大于上四分位数和小于下四分位数的样本
    selected_samples_lower <- as.matrix(names(df_gene)[df_gene < lower_quartile])
    rownames(selected_samples_lower) = selected_samples_lower
    selected_samples_lower[,] = rep(1, length(selected_samples_lower))
    selected_samples_lower = as.matrix(apply(selected_samples_lower,1,as.integer))
    
    selected_samples_up <- as.matrix(names(df_gene)[df_gene > upper_quartile ])
    rownames(selected_samples_up) = selected_samples_up
    selected_samples_up[,] = rep(2, length(selected_samples_up))
    selected_samples_up = as.matrix(apply(selected_samples_up,1,as.integer))
    if (length(selected_samples_lower) == 0){
      selected_samples_lower = selected_samples_up
      selected_samples_lower[,1] = 1
      print("#######################selected_samples_lower is error################################")
      print(rownames(df_gene))
    }
    if (length(selected_samples_up) == 0){
      selected_samples_up = selected_samples_lower
      selected_samples_up[,1] = 2
      print("#######################selected_samples_up is error################################")
      print(rownames(df_gene))
    }
    cluster_output = rbind(selected_samples_lower,selected_samples_up)
    survival_result[[gene]]=SurvivalPack(Survival_Information=Survival_Information,ClusterPackResult=cluster_output)
  }
  
  return(survival_result)
  
}
#生存单基因共表达算法
Survival_modue = function(Data,Survival_Information,node_col){
  library(readxl)
  Data = Data
  Survival_Information = Survival_Information
  node_col = node_col
  df1 = as.data.frame(Data)
  df2 = as.data.frame(node_col)
  df <- df1[rownames(df1) %in% df2$name,]
  #1.按照单基因排序，筛选出大于上四分位数和小于下四分位数的样本
  survival_result_modue = NULL
  f<-function(vec,i)
  {
    X <- vec[i]
    Xn <- vec[-i]
    sum((X - Xn) ^ 2)
  }
  for (gene in 1:nrow(df)){
    df_gene <-df[gene,]
    #modue算法
    distance <-apply(df, 2, f, i=gene)
    upper_quartile_modue <- quantile(distance, 0.75)
    lower_quartile_modue <- quantile(distance, 0.25)
    selected_samples_up_modue <- as.matrix(names(distance)[distance > upper_quartile_modue])
    rownames(selected_samples_up_modue) = selected_samples_up_modue
    selected_samples_up_modue[,] = rep(1, length(selected_samples_up_modue))
    selected_samples_up_modue = as.matrix(apply(selected_samples_up_modue,1,as.integer))
    
    selected_samples_lower_modue <- as.matrix(names(distance)[distance < lower_quartile_modue])
    rownames(selected_samples_lower_modue) = selected_samples_lower_modue
    selected_samples_lower_modue[,] = rep(2, length(selected_samples_lower_modue))
    selected_samples_lower_modue = as.matrix(apply(selected_samples_lower_modue,1,as.integer))
    if (length(selected_samples_lower_modue) == 0){
      selected_samples_lower_modue = selected_samples_up_modue
      selected_samples_lower_modue[,1] = 1
      print("#######################selected_samples_lower_modue is error################################")
      print(rownames(df_gene))
    }
    if (length(selected_samples_up_modue) == 0){
      selected_samples_up_modue = selected_samples_lower_modue
      selected_samples_up_modue[,1] = 2
      print("#######################selected_samples_up_modue is error################################")
      print(rownames(df_gene))
    }
    cluster_output_modue = rbind(selected_samples_up_modue, selected_samples_lower_modue)
    survival_result_modue[[gene]]=SurvivalPack(Survival_Information=Survival_Information,ClusterPackResult=cluster_output_modue)
  }
  return(survival_result_modue)
}
#获得基因名称
Survival_gene_name = function(Data,node_col){
  Data = Data
  node_col = node_col
  df1 = as.data.frame(Data)
  df2 = as.data.frame(node_col)
  df <- df1[rownames(df1) %in% df2$name,]
  gene_name = NULL
  for (gene in 1:nrow(df)){
    df_gene <-df[gene,]
    # 获得基因的名称
    gene_name[[gene]] = rownames(df_gene)
  }
  return(gene_name)
}

# dir.create("survival")
# for (g in 1:length(Survival_gene_name)){
#   dir.create(paste("survival/",Survival_gene_name[[g]],sep=""))
#   survival_result_gene = survival_result[[g]]
#   survival_result_modue_gene = survival_result_modue[[g]]
#   write.table(survival_result_gene$overall_pvalue,file = paste("survival/",Survival_gene_name[[g]],"/overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
#   write.table(survival_result_modue_gene$overall_pvalue,file = paste("survival/",Survival_gene_name[[g]],"/modue_overall_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F,col.names = F)
#   write.table(survival_result_gene$pairwise_pvalue,file = paste("survival/",Survival_gene_name[[g]],"/Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
#   write.table(survival_result_modue_gene$pairwise_pvalue,file = paste("survival/",Survival_gene_name[[g]],"/modue_Pvalue.txt",sep=""),quote = F,sep = "\t",row.names = F)
#   for (j in 1:length(survival_result_gene$Survival_Information_json_list))
#   {
#     write.table(survival_result_gene$Survival_Information_json_list[[j]],file = paste("survival/",Survival_gene_name[[g]],"/survival_",j,".txt",sep=""),quote = F,sep = "\t",row.names = F)
#     write.table(survival_result_modue_gene$Survival_Information_json_list[[j]],file = paste("survival/",Survival_gene_name[[g]],"/modue_survival_",j,".txt",sep=""),quote = F,sep = "\t",row.names = F)
#   }
# }



