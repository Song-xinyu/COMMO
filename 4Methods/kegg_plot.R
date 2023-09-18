# KEGG富集分析气泡图
# data:KEGG富集分析结果(csv/txt/xlsx)
# 必须包含名为Description/p.adjust/Count/GeneRatio的四列
# 绘制Top10/20气泡图
# filename:输出文件名
# dpi:设置图片分辨率
# format:设置图片输出格式(png/jpg/jpeg/tiff/pdf)

plot_keggdot = function(data,num = 10,filename = "kegg_dotplot",
                        dpi = 600,format = "png",type = "gene"){
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(stringr)
  library(dplyr)
  library(readxl)
  
  # 获取KEGG富集分析结果
  # if(grepl(".txt",data)){
  #   data = read.table(data,header = T,sep = "\t",quote = "",stringsAsFactors = F)
  # }else if(grepl(".csv",data)){
  #   data = read.csv(data,stringsAsFactors = F)
  # }else if(grepl(".xlsx",data)){
  #   data = data.frame(read_xlsx(data))
  # }
  
  # 输出结果Top10/20/30绘制气泡图(默认按p.adjust升序排序)
  if(nrow(data) > 0){
    if(num > nrow(data)){
      kegg.df = data
    }else{
      kegg.df = data[1:num,]
    }
    
    # 设置图片大小
    if(nrow(kegg.df) <= 10){
      h = 5
    }else if(nrow(kegg.df) > 10 & nrow(kegg.df) <= 20){
      h = 7
    }else if(nrow(kegg.df) > 30 & nrow(kegg.df) <= 30){
      h = 9
    }
    
    # 根据type设置绘图
    if(type == "gene"){
      if(class(kegg.df$GeneRatio) == "numeric"){
        ratio = kegg.df$GeneRatio
      }else{
        gene.num = as.numeric(unlist(strsplit(kegg.df$GeneRatio[1],split = "/"))[2]) # 获取差异基因数目
        ratio = kegg.df$Count/gene.num # 计算generatio用于绘图横坐标
      }
      kegg.df$ratio = ratio
      x.label = "GeneRatio" # x轴名称
    }else if(type == "meta"){
      if(class(kegg.df$MetabolismRatio) == "numeric"){
        ratio = kegg.df$MetabolismRatio
      }else{
        metabolism.num = as.numeric(unlist(strsplit(kegg.df$MetabolismRatio[1],split = "/"))[2]) # 获取差异代谢物数目
        ratio = kegg.df$Count/metabolism.num # 计算metabolismratio用于绘图横坐标
      }
      kegg.df$ratio = ratio
      x.label = "MetabolismRatio" # x轴名称
    }
    
    # 使画出的KEGG term与输入顺序一致
    kegg.df$Description = factor(kegg.df$Description,levels = rev(kegg.df$Description))
    
    # reorder使纵轴按照term和count排序
    kegg.dot = ggplot(data = kegg.df,aes(x = ratio, y = reorder(Description,Count)))+
      geom_point(aes(size = Count,color = -log10(p.adjust)))+
      theme_bw()+
      scale_colour_gradient(low = "purple",high = "red")+
      scale_y_discrete(labels = function(x) str_wrap(x,width = 40))+
      labs(x = "GeneRatio",y = "KEGG pathway",title = "Dotplot of KEGG Pathways",
           color = expression(-log10(p.adjust)),size = "Count")+
      theme(axis.title = element_text(size = 13),axis.text = element_text(size = 11),
            plot.title = element_text(size = 15,hjust = 0.5,face = "bold"),
            legend.title = element_text(size = 12),legend.text = element_text(size = 10),
            plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
    
    # 返回结果
    return(kegg.dot)
    #ggsave(kegg.dot,filename = paste(filename,format,sep = "."),dpi = dpi,height = h)
  }else{
    message("--> plot_keggdot: No KEGG enrichment analysis results....")
  }
}