# GO富集分析条形图
# BP/CC/MF各选前十绘制(结果默认按照p.adjust升序排序)

# data: GO富集分析结果(csv/txt/xlsx)
#       必须包含名为Description/p.adjust/Count/GeneRatio/ONTOLOGY的五列, 且对应的值正确
# ont:指定绘制BP/CC/MF/ALL条形图
# num:指定用于绘制条形图的GO term数量(10/20/30)
# filename:指定输出文件名
# dpi:指定图片分辨率(300/600)
# format:指定输出图片格式(png/jpg/jpeg/tiff/pdf/svg)

gobarplot = function(data,ont = "ALL",num = 10,
                     filename = "go_barplot",dpi = 600,format = "png"){
  library(readxl)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(stringr)
  
  # 读取差异分析结果
  # if(grepl(".txt",data)){
  #   data <- read.table(data,header = T,sep = "\t",quote = "",stringsAsFactors = F)
  # }else if(grepl(".csv",data)){
  #   data <- read.csv(data,stringsAsFactors = F)
  # }else if(grepl(".xlsx",data)){
  #   data <- data.frame(read_xlsx(data))
  # }
  
  if(nrow(data) > 0){
    if(ont == "BP"){
      if(sum(data$ONTOLOGY == "BP") > 0){
        if(sum(data$ONTOLOGY == "BP") > num){
          go.df = data[data$ONTOLOGY == "BP",][1:num,]
        }else{
          go.df = data[data$ONTOLOGY == "BP",]
        }
      }else{
        go.df = data.frame()
      }
    }else if(ont == "CC"){
      if(sum(data$ONTOLOGY == "CC") > 0){
        if(sum(data$ONTOLOGY == "CC") > num){
          go.df = data[data$ONTOLOGY == "CC",][1:num,]
        }else{
          go.df = data[data$ONTOLOGY == "CC",]
        }
      }else{
        go.df = data.frame()
      }
    }else if(ont == "MF"){
      if(sum(data$ONTOLOGY == "MF") > 0){
        if(sum(data$ONTOLOGY == "MF") > num){
          go.df = data[data$ONTOLOGY == "MF",][1:num,]
        }else{
          go.df = data[data$ONTOLOGY == "MF",]
        }
      }else{
        go.df = data.frame()
      }
    }else if(ont == "ALL"){
      # BP
      if(sum(data$ONTOLOGY == "BP") > 0){
        if(sum(data$ONTOLOGY == "BP") > 10){
          go.bp = data[data$ONTOLOGY == "BP",][1:10,]
        }else{
          go.bp = data[data$ONTOLOGY == "BP",]
        }
      }else{
        go.bp = data.frame()
      }
      # CC
      if(sum(data$ONTOLOGY == "CC") > 0){
        if(sum(data$ONTOLOGY == "CC") > 10){
          go.cc = data[data$ONTOLOGY == "CC",][1:10,]
        }else{
          go.cc = data[data$ONTOLOGY == "CC",]
        }
      }else{
        go.cc = data.frame()
      }
      # MF
      if(sum(data$ONTOLOGY == "MF") > 0){
        if(sum(data$ONTOLOGY == "MF") > 10){
          go.mf = data[data$ONTOLOGY == "MF",][1:10,]
        }else{
          go.mf = data[data$ONTOLOGY == "MF",]
        }
      }else{
        go.mf = data.frame()
      }
      # 合并数据框
      go.df = rbind(go.bp,go.cc,go.mf)
    }
    
    if(nrow(go.df) > 0){
      # 设置图片尺寸
      if(nrow(go.df) <= 10){
        h = 5
        w = 7
      }else if(nrow(go.df) > 10 & nrow(go.df) <= 15){
        h = 6
        w = 7
      }else if(nrow(go.df) > 15 & nrow(go.df) <= 20){
        h = 7
        w = 8
      }else if(nrow(go.df) > 20 & nrow(go.df) <= 25){
        h = 8
        w = 9
      }else if(nrow(go.df) > 25 & nrow(go.df) <= 30){
        h = 10
        w = 9
      }
      
      # 使绘图条目顺序一致
      go.df$Description = factor(go.df$Description,levels = rev(go.df$Description))
      
      # 绘图
      p.bar = ggplot(data = go.df,aes(x = Description,y = Count,fill = ONTOLOGY))+
        geom_bar(stat = "identity",width = 0.85)+
        coord_flip()+theme_bw()+
        scale_x_discrete(labels = function(x) str_wrap(x,width = 55))+
        labs(x = "GO terms",y = "GeneNumber",title = "Barplot of GO terms")+
        theme(axis.title = element_text(size = 12),axis.text = element_text(size = 10),
              plot.title = element_text(size = 15,hjust = 0.5,face = "bold"),
              legend.title = element_text(size = 12),legend.text = element_text(size = 10),
              plot.margin = unit(c(0.4,0.4,0.4,0.4),"cm"))
      
      # 输出结果
      #ggsave(p.bar,filename = paste(filename,format,sep = "."),dpi = dpi,width = w,height = h)
      return(p.bar)
    }else{
      message("--> gobarplot: No go.df....")
    }
  }else{
    message("--> gobarplot: No GO enrichment analysis results....")
  }
}


