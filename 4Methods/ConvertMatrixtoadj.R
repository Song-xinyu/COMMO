##20171129byhes，本函数是将邻接矩阵转化为三元邻接表，第一列为列名，第二列为行名，第三列为值（包括0）
ConvertMatrixtoadj<- function(ADMatrix)
{
  temp=cbind(rep(colnames(ADMatrix),each=nrow(ADMatrix)),rep(rownames(ADMatrix),ncol(ADMatrix)))
  result=data.frame(temp,as.vector(ADMatrix))
  colnames(result)=c("colname","rowname","value")
  return(result)
}


