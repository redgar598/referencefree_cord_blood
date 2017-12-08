####################################################################
## post hoc smoking addition
####################################################################


## load in testing data
load("~/ewas3rdround/WB_betas_BMIQ_comabt_alloutliersremoved.rdata")
##meg has kindly already filtered the WB object
beta_variable<-as.data.frame(validation_betas.combat)
dim(beta_variable)

library(ggplot2)
library(reshape)
library(RCurl)

## need mval for reference free
Mval<-function(beta) log2(beta/(1-beta))

mval_variable = apply(as.data.frame(beta_variable), 1, Mval) # need mvalues for combat
mval_variable = as.data.frame(mval_variable)
mval_variable = t(mval_variable)

## impute missing mval with row median
NA2med <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
beta_variable<-t(apply(beta_variable,1, NA2med))
mval_variable<-t(apply(mval_variable,1, NA2med))


## load meta data
meta_cord<-read.table("~/ewas3rdround/AnalysisfileUBC.dat", sep="\t", header=T)
meta_cord<-meta_cord[which(meta_cord$Sample_ID%in%colnames(beta_variable)),]
meta_cord<-meta_cord[match(colnames(beta_variable), meta_cord$Sample_ID),]
identical(colnames(beta_variable), as.character(meta_cord$Sample_ID))

meta_cord$GENDER<-as.factor(meta_cord$GENDER)




################### SVA Smoke
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GENDER)
mod0 = model.matrix(~1, data.frame(meta_cord$GENDER))

## surrogates
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
sv_unsup_smoke<-svobj$sv

################### SVA supervised Smoke
library(sva)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%probes),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GENDER)
mod0 = model.matrix(~1, data.frame(meta_cord$GENDER))

## surrogates
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_smoke<-svobj$sv

################### RUV Smoke
library(missMethyl)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")

# which probes define cell type
ctl<-(rownames(mval_variable)%in%probes)

fit = RUVfit(data=mval_variable, design=as.numeric(meta_cord$GENDER),  ctl=ctl, k=5, method=c("ruv4"))
ruv_smoke<-t(fit$W)

save(sv_unsup_smoke, sv_sup_smoke, ruv_smoke, file="~/ewas3rdround/Components_smoke.Rdata")





#####################
load("~/ewas3rdround/Components_smoke.Rdata")


##################### Corrected Betas SVA unsup smoke
sv_unsup_smoke<-as.data.frame(sv_unsup_smoke)
rownames(sv_unsup_smoke)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_unsup_smoke) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.unsup.smoke <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.unsup.smoke,file="~/ewas3rdround/adj.residuals_sva.unsup.smoke.Rdata")

rm(adj.residuals.sva.unsup.smoke)
rm(betas.lm)
gc()

##################### Corrected Betas SVA sup smoke
sv_sup_smoke<-as.data.frame(sv_sup_smoke)
rownames(sv_sup_smoke)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_sup_smoke) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.sup.smoke <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.sup.smoke,file="~/ewas3rdround/adj.residuals_sva.sup.smoke.Rdata")

rm(adj.residuals.sva.sup.smoke)
rm(betas.lm)
gc()

##################### Corrected Betas RUV smoke
ruv_smoke<-as.data.frame(ruv_smoke)
rownames(ruv_smoke)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=ruv_smoke) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.ruv.smoke <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.smoke,file="~/ewas3rdround/adj.residuals_ruv.smoke.Rdata")

rm(adj.residuals.ruv.smoke)
rm(betas.lm)
gc()