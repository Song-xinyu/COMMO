##########################################################test
# setwd("D:/job/jobs/comsuc2.0/run_comsuc2/")
## TCGA
# cancer_type="BLCA"
# projects="TCGA"
## CPTAC
# projects="CPTAC"
# CPTAC_type="breast"
# omics_wd="D:/job/jobs/comsuc2.0/run_comsuc2/userdata/User_mRNA.txt"
# omics_wd="D:/job/jobs/comsuc2.0/run_comsuc2/userdata/User_CPTAC.txt"
#########################################################

###加载包

setwd("/comsuc/run_comsuc2")
source("./4Methods/4.0FunctionPack.R")
source("./3FeatureSelection/3.0FeratureSelectionFunction.R")
source("./2Platform/2.1Batch_Correction.R")

#输入参数
argv=commandArgs(TRUE)
cancer_type=as.character(argv[1]) #which type of cancer?
projects=as.character(argv[2])#which project data?
omics_wd=NULL
omics_wd=argv[3]#user file path? 

#如果用户没有上传文件则默认服务器的备用数据
if(mode(omics_wd) == "NULL"){
  if (projects=="TCGA"){
    omics_wd="D:/job/jobs/comsuc2.0/run_comsuc2/userdata/User_mRNA.txt"
  }
  if (projects=="CPTAC"){
    omics_wd="D:/job/jobs/comsuc2.0/run_comsuc2/userdata/User_CPTAC.txt"
  }
}

#加载本地数据库
if (projects=="TCGA")
{
  omics=c(1,0)
  names(omics)=c("TCGA","CPTAC")
  TCGA_sta=as.matrix(read.table(file = "../run_comsuc2/publicData/TCGA_sta.txt",header=T,sep="\t",quote="",check.names = F))
  cancer_index=which(rownames(TCGA_sta)==cancer_type)
  
  TCGA_list=NULL
  if (TCGA_sta[cancer_index,4]!=0)#mRNA
  {
    TCGA_list[[1]]=as.matrix(read.table(file=paste("../run_comsuc2/publicData/TCGA_preprocess1.1/TCGA_UCSC/",cancer_index,cancer_type,"/HiSeqV2_PANCAN",sep=""),header=T,sep="\t",quote="",check.names = F))
  }
}

if (projects=="CPTAC")
{
  omics=c(0,1)
  names(omics)=c("TCGA","CPTAC")
  CPTAC_list=NULL
  if (CPTAC_type!="N")
  {
    CPTAC_list[[1]]=as.matrix(read.table(file=paste("./publicData/CPTAC/",CPTAC_type,"/type",sep=""),header=T,sep="\t",quote="",check.names = F))
  }
}

#加载User数据
User_list=NULL
if (omics[1]==1)#mRNA
{
  User_list[[1]]=as.matrix(read.table(file = omics_wd, header=T,sep="\t",quote="",check.names = F))
}
if (omics[2]==1)#CPTAC
{
  User_list[[2]]=as.matrix(read.table(file = omics_wd, header=T,sep="\t",quote="",check.names = F))
}

####3 Data preprocessing####
Background_list=NULL
if (projects=="TCGA")
{
  Background_list[[1]]=TCGA_list[[1]]
}
if (projects=="CPTAC")
{
  Background_list[[2]]=CPTAC_list[[1]]
}


for (i in 1:length(omics))#留下用户和本地数据库交集的数据
{
  if (omics[i]!=0)
  {
    name1=rownames(Background_list[[i]])
    name2=rownames(User_list[[i]])
    name_overlap=intersect(name1,name2)
    Background_list[[i]]=Background_list[[i]][name_overlap,]
    User_list[[i]]=User_list[[i]][name_overlap,]
    rm(name1,name2,name_overlap)
  }
}
####3.2 Remove batch effect, combine project data and user data####
Data_List=NULL

if (omics[1]!=0)
{
  Data_List[[1]]=Batch_Correction(Background_list[[1]],User_list[[1]])
}
if (omics[2]!=0)
{
  Data_List[[1]]=Batch_Correction(Background_list[[2]],User_list[[2]])
}

#输出聚类算法的输入文件
#dim(Data_List[[1]])
Data_shape=dim(Data_List[[1]])[2]
write.table(Data_shape,file = "./User_datashape.txt",quote = F,row.names = F,col.names = F)
write.table(Data_List[[1]], file = "./Input_file.txt", quote = F,sep = "\t",col.names = T,row.names = T)

