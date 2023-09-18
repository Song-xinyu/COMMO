####1 Set argument, source RScript, and library packages####
argv=commandArgs(TRUE)  
user_files=argv[1]
#test data
setwd("D:/job/jobs/comsuc2.0/run_comsuc2/")
user_files="D:/job/jobs/comsuc2.0/run_comsuc2/publicData/CPTAC/breast/type"

####2 Load Input Data####
User_list=NULL
User_list[[1]]=as.matrix(read.table(file = user_files ,header=T,sep="\t",quote="",check.names = F))

####3 Data Check Function####
omics_index=1
DataCheck=function(User_list){
  message_txt=NULL
  for (i in omics_index)
  {
    if (mode(User_list[[i]])!="numeric")
    {
      message_txt="Error! There are some non-numeric value in omics data. Please check and re-upload."
      return(message_txt)
    }
  }
  for (i in omics_index)
  {
    if (length(which(is.na(User_list[[i]])))!=0)
    {
      message_txt="Error! There are some NULL value in omics data. Please check and re-upload."
      return(message_txt)
    }
  }

  return(message_txt)
}

####4 Data Check####
user_file = basename(user_files)
message_output=DataCheck(User_list)

if (is.null(message_output))
{

  Data_shape=NULL
  Data_shape=length(colnames(User_list[[1]]))
  write.table(Data_shape,file = "./User_datashape.txt",quote = F,row.names = F,col.names = F)
  write.table(User_list[[1]],file = "./User_data.txt",quote = F,row.names = T,col.names = T,sep="\t")

}else{
  write.table(message_output,file = "./Upload_error.txt",quote = F,row.names = F,col.names = F)
}

