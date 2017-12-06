###############################################################################################
##MAE
###############################################################################################
library(hydroGOF)


MSE_compare_all<-function(betas, method2_name){

    ## function to get MAE for each dataframe
    method_MAE<-function(adj.res, method_name){
      metod_beta<-as.data.frame(adj.res)
      cpg_mae<-mae(t(metod_beta),t(betas))
      data.frame(method_1=method_name,method_2=method2_name, mn_cpg_mae=mean(cpg_mae))
    }
    
    ## load GS
    load("~/ewas3rdround/facs_pcs_corrected_betas.Rdata")
    gold_standard<-as.data.frame(adj.residuals)
    error_GS<-method_MAE(gold_standard, "FACS - PCA - Gold-Standard")
    
    ## Decon PCA
    load("~/ewas3rdround/decon_pcs_corrected_betas.Rdata")
    err<-method_MAE(adj.residuals, "Deconvolution - PCA")
    error_GS<-rbind(error_GS, err)
    
    ## Decon counts
    load("~/ewas3rdround/decon_corrected_betas.Rdata")
    err<-method_MAE(adj.residuals, "Deconvolution - Drop One Cell Type")
    error_GS<-rbind(error_GS, err)
    
    ## FACS PCA
    load("~/ewas3rdround/facs_corrected_betas.Rdata")
    err<-method_MAE(adj.residuals, "FACS - Drop One Cell Type")
    error_GS<-rbind(error_GS, err)
    
    
    ## Uncorrected
    load("~/ewas3rdround/WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
    err<-method_MAE(validation_betas.combat, "Uncorrected")
    error_GS<-rbind(error_GS, err)
    
    ## Refactor
    load("~/ewas3rdround/adj.residuals_refactor.Rdata")
    err<-method_MAE(adj.residuals.refactor, "ReFACTor")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.refactor)
    
    ## reffreecellmix
    load("~/ewas3rdround/adj.residuals_reffreecellmix.Rdata")
    err<-method_MAE(adj.residuals.reffreecellmix, "RefFreeCellMix")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.reffreecellmix)
    
    ## RUV GA
    load("~/ewas3rdround/adj.residuals_ruv.ga.Rdata")
    err<-method_MAE(adj.residuals.ruv.ga, "RUV - GA")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.ruv.ga)
    
    ## RUV sex
    load("~/ewas3rdround/adj.residuals_ruv.sex.Rdata")
    err<-method_MAE(adj.residuals.ruv.sex, "RUV - Sex")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.ruv.sex)
    
    
    
    ## SVA - Supervised GA 
    load("~/ewas3rdround/adj.residuals_sva.sup.ga.Rdata")
    err<-method_MAE(adj.residuals.sva.sup.ga, "SVA - Supervised GA")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.sva.sup.ga)
    
    ## SVA - Unsupervised GA 
    load("~/ewas3rdround/adj.residuals_sva.unsup.ga.Rdata")
    err<-method_MAE(adj.residuals.sva.unsup.ga, "SVA - Unsupervised GA")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.sva.unsup.ga)
    
    ## SVA - Supervised Sex 
    load("~/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")
    err<-method_MAE(adj.residuals.sva.sup.sex, "SVA - Supervised Sex")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.sva.sup.sex)
    
    ## SVA - Unsupervised Sex 
    load("~/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")
    err<-method_MAE(adj.residuals.sva.unsup.sex, "SVA - Unsupervised Sex")
    error_GS<-rbind(error_GS, err)
    rm(adj.residuals.sva.unsup.sex)
  
    error_GS}

# load("~/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")
# method_adj<-as.data.frame(adj.residuals.sva.unsup.sex)
# MSE_compare_all(method_adj, "SVA - Unsupervised Sex")


## GS
load("~/ewas3rdround/facs_pcs_corrected_betas.Rdata")
method_adj<-as.data.frame(adj.residuals)
error_all<-MSE_compare_all(method_adj, "FACS - PCA - Gold-Standard")

## Decon PCA
load("~/ewas3rdround/decon_pcs_corrected_betas.Rdata")
method_adj<-as.data.frame(adj.residuals)
err<-MSE_compare_all(method_adj, "Deconvolution - PCA")
error_all<-rbind(error_all, err)

## Decon counts
load("~/ewas3rdround/decon_corrected_betas.Rdata")
method_adj<-as.data.frame(adj.residuals)
err<-MSE_compare_all(method_adj, "Deconvolution - Drop One Cell Type")
error_all<-rbind(error_all, err)

## FACS PCA
load("~/ewas3rdround/facs_corrected_betas.Rdata")
method_adj<-as.data.frame(adj.residuals)
err<-MSE_compare_all(method_adj, "FACS - Drop One Cell Type")
error_all<-rbind(error_all, err)

## Uncorrected
load("~/ewas3rdround/WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
method_adj<-as.data.frame(validation_betas.combat)
err<-MSE_compare_all(method_adj, "Uncorrected")
error_all<-rbind(error_all, err)

## Refactor
load("~/ewas3rdround/adj.residuals_refactor.Rdata")
method_adj<-as.data.frame(adj.residuals.refactor)
err<-MSE_compare_all(method_adj, "ReFACTor")
error_all<-rbind(error_all, err)

## reffreecellmix
load("~/ewas3rdround/adj.residuals_reffreecellmix.Rdata")
method_adj<-as.data.frame(adj.residuals.reffreecellmix)
err<-MSE_compare_all(method_adj, "RefFreeCellMix")
error_all<-rbind(error_all, err)

## RUV GA
load("~/ewas3rdround/adj.residuals_ruv.ga.Rdata")
method_adj<-as.data.frame(adj.residuals.ruv.ga)
err<-MSE_compare_all(method_adj, "RUV - GA")
error_all<-rbind(error_all, err)

## RUV sex
load("~/ewas3rdround/adj.residuals_ruv.sex.Rdata")
method_adj<-as.data.frame(adj.residuals.ruv.sex)
err<-MSE_compare_all(method_adj, "RUV - Sex")
error_all<-rbind(error_all, err)

## SVA - Supervised GA 
load("~/ewas3rdround/adj.residuals_sva.sup.ga.Rdata")
method_adj<-as.data.frame(adj.residuals.sva.sup.ga)
err<-MSE_compare_all(method_adj, "SVA - Supervised GA")
error_all<-rbind(error_all, err)

## SVA - Unsupervised GA 
load("~/ewas3rdround/adj.residuals_sva.unsup.ga.Rdata")
method_adj<-as.data.frame(adj.residuals.sva.unsup.ga)
err<-MSE_compare_all(method_adj, "SVA - Unsupervised GA")
error_all<-rbind(error_all, err)

## SVA - Supervised Sex 
load("~/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")
method_adj<-as.data.frame(adj.residuals.sva.sup.sex)
err<-MSE_compare_all(method_adj, "SVA - Supervised Sex")
error_all<-rbind(error_all, err)

## SVA - Unsupervised Sex 
load("~/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")
method_adj<-as.data.frame(adj.residuals.sva.unsup.sex)
err<-MSE_compare_all(method_adj, "SVA - Unsupervised Sex")
error_all<-rbind(error_all, err)


save(error_all, file="~/RE_GenR/error_all.RData")









###############################################################################################
## Correlation
###############################################################################################
## comapre all to gold standard
cor_GS<-function(betas, method_name){
  cpg_cor<-sapply(1:nrow(gold_standard), function(x){
  cor(as.numeric(gold_standard[x,]),betas[x,])})
  data.frame(correlation=cpg_cor, method=method_name)}

## load GS
load("~/ewas3rdround/facs_pcs_corrected_betas.Rdata")
gold_standard<-as.data.frame(adj.residuals)

## Decon PCA
load("~/ewas3rdround/decon_pcs_corrected_betas.Rdata")
cor_to_GS<-cor_GS(adj.residuals, "Deconvolution - PCA")

## Decon counts
load("~/ewas3rdround/decon_corrected_betas.Rdata")
err<-cor_GS(adj.residuals, "Deconvolution - Drop One Cell Type")
cor_to_GS<-rbind(cor_to_GS, err)

## FACS PCA
load("~/ewas3rdround/facs_corrected_betas.Rdata")
err<-cor_GS(adj.residuals, "FACS - Drop One Cell Type")
cor_to_GS<-rbind(cor_to_GS, err)

## Uncorrected
load("~/ewas3rdround/WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
err<-cor_GS(validation_betas.combat, "Uncorrected")
cor_to_GS<-rbind(cor_to_GS, err)

## Refactor
load("~/ewas3rdround/adj.residuals_refactor.Rdata")
err<-cor_GS(adj.residuals.refactor, "ReFACTor")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.refactor)

## reffreecellmix
load("~/ewas3rdround/adj.residuals_reffreecellmix.Rdata")
err<-cor_GS(adj.residuals.reffreecellmix, "RefFreeCellMix")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.reffreecellmix)

## RUV GA
load("~/ewas3rdround/adj.residuals_ruv.ga.Rdata")
err<-cor_GS(adj.residuals.ruv.ga, "RUV - GA")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.ruv.ga)

## RUV sex
load("~/ewas3rdround/adj.residuals_ruv.sex.Rdata")
err<-cor_GS(adj.residuals.ruv.sex, "RUV - Sex")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.ruv.sex)

## SVA - Supervised GA 
load("~/ewas3rdround/adj.residuals_sva.sup.ga.Rdata")
err<-cor_GS(adj.residuals.sva.sup.ga, "SVA - Supervised GA")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.sva.sup.ga)

## SVA - Unsupervised GA 
load("~/ewas3rdround/adj.residuals_sva.unsup.ga.Rdata")
err<-cor_GS(adj.residuals.sva.unsup.ga, "SVA - Unsupervised GA")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.sva.unsup.ga)

## SVA - Supervised Sex 
load("~/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")
err<-cor_GS(adj.residuals.sva.sup.sex, "SVA - Supervised Sex")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.sva.sup.sex)

## SVA - Unsupervised Sex 
load("~/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")
err<-cor_GS(adj.residuals.sva.unsup.sex, "SVA - Unsupervised Sex")
cor_to_GS<-rbind(cor_to_GS, err)
rm(adj.residuals.sva.unsup.sex)

save(cor_to_GS, file="~/RE_GenR/correlation_distributions.RData")

