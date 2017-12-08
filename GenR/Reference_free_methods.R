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


              
################# RefFreeCellMix
library(RefFreeEWAS)
#Kˆ=3K^=3 for the cord blood data set BL-as, larger for the other three peripheral blood datasets
# We performed reference-free deconvolution with the method RefFreeCellMix by Houseman et al. [42] 
# using the R package RefFreeEWAS. In accordance with the original publication of the method [42], 
# we applied it to the 20,000 most variable CpG positions from the methylation matrix, unless the 
# total number of rows was less, in which case we used the full matrix. In the former case, we used 
# the available option to obtain the estimates of the methylation components for all CpGs as the 
# final step of the deconvolution procedure (supplying the complete data matrix as argument Yfinal).

testArray1 <- RefFreeCellMix(as.matrix(beta_variable),K=5,iters=5,Yfinal=NULL)#verbose=T
RefFreeCounts<-as.data.frame(testArray1$Omega)


################### ReFACTor
#https://github.com/cozygene/refactor/tree/master/R
source("~/ewas3rdround/refactor.R")

  # only need to run this one so commented out for debugging
  beta_df<-as.data.frame(beta_variable)
  beta_df$ID<-rownames(beta_df)
  beta_df<-beta_df[,c(ncol(beta_df), 1:(ncol(beta_df)-1))]
  write.table(beta_df, file="~/ewas3rdround/betas_for_refactor.txt", quote=F, row.names=F,sep="\t")

k = 5
datafile = "~/ewas3rdround/betas_for_refactor.txt"

results <- refactor(datafile,k)
RC <- as.data.frame(results$refactor_components) # Extract the ReFACTor components


## save for plotting
save(RefFreeCounts, RC, file="~/ewas3rdround/Count_like_data.Rdata")



################### SVA GA
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GESTBIR)
mod0 = model.matrix(~1, data.frame(meta_cord$GESTBIR))

## surrogatesref
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
sv_unsup_gestage<-svobj$sv

################### SVA supervised GA
library(sva)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%probes),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GESTBIR)
mod0 = model.matrix(~1, data.frame(meta_cord$GESTBIR))

## surrogates
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_gestage<-svobj$sv

################### RUV GA
library(missMethyl)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")

# which probes define cell type
ctl<-(rownames(mval_variable)%in%probes)

fit = RUVfit(data=mval_variable, design=meta_cord$GESTBIR,  ctl=ctl, k=5, method=c("ruv4"))
ruv_GA<-t(fit$W)

save(sv_unsup_gestage, sv_sup_gestage, ruv_GA, file="~/ewas3rdround/Components_GA.Rdata")





################### SVA Sex
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GENDER)
mod0 = model.matrix(~1, data.frame(meta_cord$GENDER))

## surrogates
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
sv_unsup_sex<-svobj$sv

################### SVA supervised Sex
library(sva)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%probes),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GENDER)
mod0 = model.matrix(~1, data.frame(meta_cord$GENDER))

## surrogates
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_sex<-svobj$sv

################### RUV sex
library(missMethyl)
load("~/ewas3rdround/cord_cell_type_specific_probes.rdata")

# which probes define cell type
ctl<-(rownames(mval_variable)%in%probes)

fit = RUVfit(data=mval_variable, design=as.numeric(meta_cord$GENDER),  ctl=ctl, k=5, method=c("ruv4"))
ruv_sex<-t(fit$W)

save(sv_unsup_sex, sv_sup_sex, ruv_sex, file="~/ewas3rdround/Components_sex.Rdata")




##############################################################################################################################
#### Correction
##############################################################################################################################
load("~/ewas3rdround/Count_like_data.Rdata")
################## Corrected Betas reffreecellmix
colnames(RefFreeCounts)<-paste("comp", colnames(RefFreeCounts), sep="")

betas.lm <- apply(beta_variable, 1, function(x){
  components <- RefFreeCounts[colnames(beta_variable),]
  lm(x~comp1+comp2+comp3+comp4+comp5,data=components) 
})

# extract matrix of residuals from resulting linear models
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names

# add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals.reffreecellmix <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))

# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.
# adj.residuals.reffreecellmix[adj.residuals.reffreecellmix<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
# adj.residuals.reffreecellmix[adj.residuals.reffreecellmix>1] <- 0.999 # convert any values that are greater than 1 to 0.999
save(adj.residuals.reffreecellmix,file="~/ewas3rdround/adj.residuals_reffreecellmix.Rdata")

rm(adj.residuals.reffreecellmix)
rm(betas.lm)
gc()



##################### Corrected Betas refactor
rownames(RC)<-colnames(beta_variable)

betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~PC1+PC2+PC3+PC4+PC5,data=RC) 
})

# extract matrix of residuals from resulting linear models
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names

# add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals.refactor <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))

# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.
# adj.residuals.refactor[adj.residuals.refactor<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
# adj.residuals.refactor[adj.residuals.refactor>1] <- 0.999 # convert any values that are greater than 1 to 0.999
save(adj.residuals.refactor,file="~/ewas3rdround/adj.residuals_refactor.Rdata")

rm(adj.residuals.refactor)
rm(betas.lm)
gc()



#####################
load("~/ewas3rdround/Components_GA.Rdata")


##################### Corrected Betas SVA unsup GA
sv_unsup_ga<-as.data.frame(sv_unsup_gestage)
rownames(sv_unsup_ga)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_unsup_ga) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.unsup.ga <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.unsup.ga,file="~/ewas3rdround/adj.residuals_sva.unsup.ga.Rdata")

rm(adj.residuals.sva.unsup.ga)
rm(betas.lm)
gc()

##################### Corrected Betas SVA sup GA
sv_sup_ga<-as.data.frame(sv_sup_gestage)
rownames(sv_sup_ga)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_sup_ga) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.sup.ga <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.sup.ga,file="~/ewas3rdround/adj.residuals_sva.sup.ga.Rdata")

rm(adj.residuals.sva.sup.ga)
rm(betas.lm)
gc()

##################### Corrected Betas RUV GA
ruv_GA<-as.data.frame(ruv_GA)
rownames(ruv_GA)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=ruv_GA) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.ruv.ga <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.ga,file="~/ewas3rdround/adj.residuals_ruv.ga.Rdata")

rm(adj.residuals.ruv.ga)
rm(betas.lm)
gc()



#####################
load("~/ewas3rdround/Components_sex.Rdata")


##################### Corrected Betas SVA unsup sex
sv_unsup_sex<-as.data.frame(sv_unsup_sex)
rownames(sv_unsup_sex)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_unsup_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.unsup.sex <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.unsup.sex,file="~/ewas3rdround/adj.residuals_sva.unsup.sex.Rdata")

rm(adj.residuals.sva.unsup.sex)
rm(betas.lm)
gc()

##################### Corrected Betas SVA sup sex
sv_sup_sex<-as.data.frame(sv_sup_sex)
rownames(sv_sup_sex)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_sup_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.sva.sup.sex <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.sup.sex,file="~/ewas3rdround/adj.residuals_sva.sup.sex.Rdata")

rm(adj.residuals.sva.sup.sex)
rm(betas.lm)
gc()

##################### Corrected Betas RUV sex
ruv_sex<-as.data.frame(ruv_sex)
rownames(ruv_sex)<-colnames(beta_variable)
betas.lm <- apply(beta_variable, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=ruv_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_variable) # re-name residuals columns with sample names
adj.residuals.ruv.sex <- residuals+matrix(apply(beta_variable, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.sex,file="~/ewas3rdround/adj.residuals_ruv.sex.Rdata")

rm(adj.residuals.ruv.sex)
rm(betas.lm)
gc()

