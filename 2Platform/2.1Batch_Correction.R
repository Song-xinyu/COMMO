Batch_Correction <- function(matrix1, matrix2){
  if(!require("sva")) BiocManager::install("sva")
  library(sva)
  matrix1 <- matrix1[which((apply(matrix1, 1, var)!=0) == "TRUE"), ]
  matrix2 <- matrix2[which((apply(matrix2, 1, var)!=0) == "TRUE"), ]
  features <- intersect(rownames(matrix1), rownames(matrix2))
  matrix1 <- matrix1[features, ]
  matrix2 <- matrix2[features, ]
  
  pheno <- as.data.frame(matrix(NA, nrow=ncol(matrix1)+ncol(matrix2), ncol=2))
  colnames(pheno) <- c("samples", "batch")
  pheno[,1] <- c(colnames(matrix1), colnames(matrix2))
  pheno[,2] <- c(rep(1,ncol(matrix1)), rep(2,ncol(matrix2)))
  batch <- pheno$batch
  modcombat <- model.matrix(~1, data=pheno)
  
  AllMatrix <- cbind(matrix1, matrix2)
  AllMatrix <- AllMatrix[complete.cases(AllMatrix), ]
  #using the function of Combat to remove/eliminate the batch effect
  #combat_edata = ComBat(dat=as.matrix(AllMatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F)
  combat_edata <- try(ComBat(dat=as.matrix(AllMatrix), batch=batch, mod=modcombat, par.prior=TRUE, prior.plots=F))
  if('try-error' %in% class(combat_edata))
  {
    combat_edata <- ComBat(dat=as.matrix(AllMatrix), batch=batch, mod=modcombat, par.prior=F, prior.plots=F)
  }
  return(combat_edata)
}