########################################## 
## EWAS
########################################## 

## load in testing data
load("/big_data/reffree/WB_betas_BMIQ_combat_together.rdata")
beta<-as.data.frame(validation_betas.combat)

library(ggplot2)
library(reshape)
library(RCurl)

## remove invariable probes
load("invariable_cordblood_CpGs.Rdata")
beta_variable<-beta[which(!(rownames(beta)%in%invar_in_beta_and_independent)),]#316128 

## need mval for reference free
Mval<-function(beta) log2(beta/(1-beta))
mval = apply(as.data.frame(beta), 1, Mval) # need mvalues for combat
mval = as.data.frame(mval)
mval = t(mval)

mval_variable = apply(as.data.frame(beta_variable), 1, Mval) # need mvalues for combat
mval_variable = as.data.frame(mval_variable)
mval_variable = t(mval_variable)

## impute missing mval with row median
NA2med <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
mval_complete <- t(apply(mval,1, NA2med))
beta_complete <- t(apply(beta,1, NA2med))
beta_variable<-t(apply(beta_variable,1, NA2med))
mval_variable<-t(apply(mval_variable,1, NA2med))


## load meta data
meta_cord<-read.csv("/big_data/redgar/cordblood/24_wholebloods_samplesheet.csv")
meta_cord<-meta_cord[match(colnames(beta), meta_cord$X),]
identical(colnames(beta_variable), as.character(meta_cord$X))

## Delta beta functions
delbeta_GA<-function(beta_df){
  sapply(1:nrow(beta_df), function(cpg) {
    z<-lm(unlist(beta_df[cpg,]) ~ meta_cord$GA)
    intercept=z$coefficients[1]
    slope=z$coefficients[2]
    as.numeric((slope*max(meta_cord$GA, na.rm=T)+intercept))-as.numeric((slope*min(meta_cord$GA, na.rm=T)+intercept))
})}


delbeta_sex<-function(beta_df){
  F_beta<-beta_df[,which(meta_cord$Sex=="F")]
  M_beta<-beta_df[,which(meta_cord$Sex=="M")]
  rowMeans(M_beta)-rowMeans(F_beta)
}


#########################################################   
## reference free SVA and RUV p values
#########################################################   


################### SVA GA
library(sva)
library(limma)

mod = model.matrix(~meta_cord$GA)
mod0 = model.matrix(~1, data.frame(meta_cord$GA))
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
modSv = cbind(mod,svobj$sv)
fit1 = lmFit(as.matrix(mval_variable),modSv)
fit_SVA_unsup_GA = eBayes(fit1)
fit_SVA_unsup_GA<-data.frame(CpG=rownames(topTable(fit_SVA_unsup_GA)),pval=topTable(fit_SVA_unsup_GA)[,c("P.Value")])

################### SVA supervised GA
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%diff_cpgs$CpG),]
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_gestage<-svobj$sv

################### RUV GA
library(missMethyl)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")

# which probes define cell type
ctl<-(rownames(mval_variable)%in%diff_cpgs$CpG)

fit = RUVfit(data=mval_variable, design=meta_cord$GA,  ctl=ctl, k=5, method=c("ruv4"))
ruv_GA<-t(fit$W)

save(sv_unsup_gestage, sv_sup_gestage, ruv_GA, file="/big_data/redgar/cordblood/Components_GA.Rdata")





################### SVA Sex
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$Sex)
mod0 = model.matrix(~1, data.frame(meta_cord$Sex))

## surrogates
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
sv_unsup_sex<-svobj$sv

################### SVA supervised Sex
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%diff_cpgs$CpG),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$Sex)
mod0 = model.matrix(~1, data.frame(meta_cord$Sex))

## surrogates
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_sex<-svobj$sv

################### RUV sex
library(missMethyl)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")

# which probes define cell type
ctl<-(rownames(mval_variable)%in%diff_cpgs$CpG)

fit = RUVfit(data=mval_variable, design=as.numeric(meta_cord$Sex),  ctl=ctl, k=5, method=c("ruv4"))
ruv_sex<-t(fit$W)

save(sv_unsup_sex, sv_sup_sex, ruv_sex, file="/big_data/redgar/cordblood/Components_sex.Rdata")







