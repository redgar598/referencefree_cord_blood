
##RMSE
library(hydroGOF)


## load GS
load("/big_data/reffree/facs_pcs_corrected_betas.Rdata")
gold_standard<-as.data.frame(adj.residuals)

## function to get RMSE to each dataframe
method_RMSE<-function(adj.res, method_name){
  metod_beta<-as.data.frame(adj.residuals)
  cpg_rmse<-rmse(t(gold_standard), t(metod_beta))
  data.frame(method=method_name, mn_cpg_rmse=mean(cpg_rmse))
}


## Decon PCA
load("/big_data/reffree/decon_pcs_corrected_betas.Rdata")
method_RMSE(adj.residuals, "Deconvolution - PCA")

## Decon counts
load("/big_data/reffree/decon_corrected_betas.Rdata")
method_RMSE(adj.residuals, "Deconvolution - Drop One Cell Type")

## FACS PCA
load("/big_data/reffree/facs_corrected_betas.Rdata")
method_RMSE(adj.residuals, "FACS - Drop One Cell Type")

## Uncorrected
load("/big_data/reffree/WB_betas_BMIQ_combat_together.rdata")
method_RMSE(validation_betas.combat, "Uncorrected")
