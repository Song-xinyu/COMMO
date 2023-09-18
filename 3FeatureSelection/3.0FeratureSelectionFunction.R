####Function: Pre-select features in each omics data matrix####
#Data: the input matrix(row represent feature, and column represent sample)
#method: MAD(default) or VAR
#toprank: top rank of selected features, 1500 default
Feature_selection=function(Data, method="MAD", toprank=1500)
{
  if (nrow(Data)>=toprank)
  {
    if (method=="VAR")
    {
      temp=apply(Data, 1, var)
      Data=Data[which(rank(-temp)<=toprank),]
    }else
    {
      temp=apply(Data, 1, mad)
      Data=Data[which(rank(-temp)<=toprank),]
    }
  }
  return(Data)
}

####1.1 CNV data feature selection####
CNV_feature_selection=function(Data)
{
  return(Feature_selection(Data=Data, method="VAR", toprank = 1500))
}

####1.2 Methylation data feature selection####
METH_feature_selection=function(Data)
{
  return(Feature_selection(Data=Data, method="VAR", toprank = 1500))
}

####1.3 mRNA data feature selection####
MRNA_feature_selection=function(Data)
{
  return(Feature_selection(Data=Data, method="MAD", toprank = 1500))
}

####1.4 miRNA data feature selection####
MIRNA_feature_selection=function(Data)
{
  return(Feature_selection(Data=Data, method="MAD", toprank = 647))
}

####1.5 RPPA data feature selection####
RPPA_feature_selection=function(Data)
{
  return(Data)
}






