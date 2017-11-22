              ##################### installations
              source("https://bioconductor.org/biocLite.R")
              biocLite("qvalue")
              
              install.packages("isva")
              install.packages("RefFreeEWAS")
              
              #sudo apt-get install libcurl4-openssl-dev
              #sudo apt-get install libcurl4-openssl-dev libssl-dev
              #sudo apt-get install libxml2-dev
              #sudo apt-get install libglu1-mesa-dev
              install.packages("RCurl")
              install.packages("rgl")
              install.packages("mixOmics")
              biocLite("minfi")
              #biocLite("FlowSorted.Blood.450k")
              #biocLite("IlluminaHumanMethylation450kmanifest")
              
              install.packages("lmtest")
              
              ## MeDeCom - fails because wrong c++ compilier
              #install.packages("devtools")
              #devtools:::install_github("lutsik/MeDeCom")
              
              biocLite("sva")
              
              biocLite("missMethyl")
              
              install.packages("ggsci")
              install.packages("hydroGOF")
              
              
              install.packages("ggridges")
              install.packages("kableExtra")
              install.packages("dplyr")



##########################################################################################################################

## load in testing data
load("/big_data/reffree/WB_betas_BMIQ_combat_together.rdata")
beta<-as.data.frame(validation_betas.combat)

library(ggplot2)
library(reshape)
library(RCurl)

## remove invariable probes
x <- getURL("https://raw.githubusercontent.com/redgar598/Tissue_Invariable_450K_CpGs/master/Invariant_Blood_CpGs.csv")
y <- read.csv(text = x)
beta_invariable<-beta[which(rownames(beta)%in%y$CpG),]#110228/114204 of the independent invariable sites are in beta
Variation<-function(x) {quantile(x, c(0.9), na.rm=T)[[1]]-quantile(x, c(0.1), na.rm=T)[[1]]}
beta_ref_range<-sapply(1:nrow(beta_invariable), function(x) Variation(beta_invariable[x,]))
Invariable_in_beta<-beta_invariable[which(beta_ref_range<0.05),]
invar_in_beta_and_independent<-intersect(y$CpG, rownames(Invariable_in_beta)) #106613/114204 (96.1%)
save(invar_in_beta_and_independent, file="invariable_cordblood_CpGs.Rdata")

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




              
################# RefFreeCellMix
library(RefFreeEWAS)
#Kˆ=3K^=3 for the cord blood data set BL-as, larger for the other three peripheral blood datasets
# We performed reference-free deconvolution with the method RefFreeCellMix by Houseman et al. [42] using the R package RefFreeEWAS. In accordance with the original publication of the method [42], we applied it to the 20,000 most variable CpG positions from the methylation matrix, unless the total number of rows was less, in which case we used the full matrix. In the former case, we used the available option to obtain the estimates of the methylation components for all CpGs as the final step of the deconvolution procedure (supplying the complete data matrix as argument Yfinal).

testArray1 <- RefFreeCellMix(as.matrix(beta_variable),K=5,iters=5,Yfinal=NULL)#verbose=T
RefFreeCounts<-as.data.frame(testArray1$Omega)


################### ReFACTor
#https://github.com/cozygene/refactor/tree/master/R
source("/big_data/redgar/cordblood/refactor.R")
        beta_df<-as.data.frame(beta_variable)
        beta_df$ID<-rownames(beta_df)
        beta_df<-beta_df[,c(ncol(beta_df), 1:(ncol(beta_df)-1))]
        write.table(beta_df, file="/big_data/redgar/cordblood/betas_for_refactor.txt", quote=F, row.names=F,sep="\t")

k = 5
datafile = "/big_data/redgar/cordblood/betas_for_refactor.txt"

results <- refactor(datafile,k)
RC <- as.data.frame(results$refactor_components) # Extract the ReFACTor components


## save for plotting
save(RefFreeCounts, RC, file="/big_data/redgar/cordblood/Count_like_data.Rdata")

################### SVA GA
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GA)
mod0 = model.matrix(~1, data.frame(meta_cord$GA))

## surrogates
svobj = sva(as.matrix(mval_variable),mod,mod0,n.sv=5)
sv_unsup_gestage<-svobj$sv

################### SVA supervised GA
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%diff_cpgs$CpG),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GA)
mod0 = model.matrix(~1, data.frame(meta_cord$GA))

## surrogates
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




##############################################################################################################################
#### Correction
##############################################################################################################################
load("/big_data/redgar/cordblood/Count_like_data.Rdata")
################## Corrected Betas reffreecellmix
colnames(RefFreeCounts)<-paste("comp", colnames(RefFreeCounts), sep="")

betas.lm <- apply(beta_complete, 1, function(x){
  components <- RefFreeCounts[colnames(beta_complete),]
  lm(x~comp1+comp2+comp3+comp4+comp5,data=components) 
})

# extract matrix of residuals from resulting linear models
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names

# add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals.reffreecellmix <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))

# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.
# adj.residuals.reffreecellmix[adj.residuals.reffreecellmix<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
# adj.residuals.reffreecellmix[adj.residuals.reffreecellmix>1] <- 0.999 # convert any values that are greater than 1 to 0.999
save(adj.residuals.reffreecellmix,file="/big_data/redgar/cordblood/adj.residuals_reffreecellmix.Rdata")

rm(adj.residuals.reffreecellmix)
rm(betas.lm)
gc()



##################### Corrected Betas refactor
rownames(RC)<-colnames(beta_complete)

betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~PC1+PC2+PC3+PC4+PC5,data=RC) 
})

# extract matrix of residuals from resulting linear models
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names

# add the residuals of each regression model to the mean methylation value of each probe (mean across all samples) to obtain the “adjusted” methylation data.
adj.residuals.refactor <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))

# To make sure we do not induce any NAs into the dataset when we convert the beta values back M-values (by log2 transformation), we need to ensure we do not have any corrected beta values that are greater or equal to zero or any beta values that are greater than 1.
# adj.residuals.refactor[adj.residuals.refactor<=0] <- 0.001 # convert any values that are less than or equal to zero to 0.001
# adj.residuals.refactor[adj.residuals.refactor>1] <- 0.999 # convert any values that are greater than 1 to 0.999
save(adj.residuals.refactor,file="/big_data/redgar/cordblood/adj.residuals_refactor.Rdata")

rm(adj.residuals.refactor)
rm(betas.lm)
gc()



#####################
load("/big_data/redgar/cordblood/Components_GA.Rdata")


##################### Corrected Betas SVA unsup GA
sv_unsup_ga<-as.data.frame(sv_unsup_gestage)
rownames(sv_unsup_ga)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_unsup_ga) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.sva.unsup.ga <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.unsup.ga,file="/big_data/redgar/cordblood/adj.residuals_sva.unsup.ga.Rdata")

rm(adj.residuals.sva.unsup.ga)
rm(betas.lm)
gc()

##################### Corrected Betas SVA sup GA
sv_sup_ga<-as.data.frame(sv_sup_gestage)
rownames(sv_sup_ga)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_sup_ga) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.sva.sup.ga <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.sup.ga,file="/big_data/redgar/cordblood/adj.residuals_sva.sup.ga.Rdata")

rm(adj.residuals.sva.sup.ga)
rm(betas.lm)
gc()

##################### Corrected Betas RUV GA
ruv_GA<-as.data.frame(ruv_GA)
rownames(ruv_GA)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=ruv_GA) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.ruv.ga <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.ga,file="/big_data/redgar/cordblood/adj.residuals_ruv.ga.Rdata")

rm(adj.residuals.ruv.ga)
rm(betas.lm)
gc()



#####################
load("/big_data/redgar/cordblood/Components_sex.Rdata")


##################### Corrected Betas SVA unsup sex
sv_unsup_sex<-as.data.frame(sv_unsup_sex)
rownames(sv_unsup_sex)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_unsup_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.sva.unsup.sex <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.unsup.sex,file="/big_data/redgar/cordblood/adj.residuals_sva.unsup.sex.Rdata")

rm(adj.residuals.sva.unsup.sex)
rm(betas.lm)
gc()

##################### Corrected Betas SVA sup sex
sv_sup_sex<-as.data.frame(sv_sup_sex)
rownames(sv_sup_sex)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=sv_sup_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.sva.sup.sex <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.sva.sup.sex,file="/big_data/redgar/cordblood/adj.residuals_sva.sup.sex.Rdata")

rm(adj.residuals.sva.sup.sex)
rm(betas.lm)
gc()

##################### Corrected Betas RUV sex
ruv_sex<-as.data.frame(ruv_sex)
rownames(ruv_sex)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5,data=ruv_sex) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.ruv.sex <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.sex,file="/big_data/redgar/cordblood/adj.residuals_ruv.sex.Rdata")

rm(adj.residuals.ruv.sex)
rm(betas.lm)
gc()



#########################################################################################################
## remove invariable probes
#########################################################################################################
load("invariable_cordblood_CpGs.Rdata")
beta_variable<-beta[which(!(rownames(beta)%in%invar_in_beta_and_independent)),]


## Uncorrected
load("/big_data/reffree/WB_betas_BMIQ_combat_together.rdata")
validation_betas.combat<-validation_betas.combat[which(!(rownames(validation_betas.combat)%in%invar_in_beta_and_independent)),] 
save(validation_betas.combat, file="/big_data/reffree/WB_betas_BMIQ_combat_together_invar_filtered.rdata")

## Refactor
load("/big_data/reffree/adj.residuals_refactor.Rdata")
adj.residuals.refactor<-adj.residuals.refactor[which(!(rownames(adj.residuals.refactor)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.refactor, file="/big_data/reffree/adj.residuals_refactor_invar_filtered.Rdata")

## reffreecellmix
load("/big_data/reffree/adj.residuals_reffreecellmix.Rdata")
adj.residuals.reffreecellmix<-adj.residuals.reffreecellmix[which(!(rownames(adj.residuals.reffreecellmix)%in%invar_in_beta_and_independent)),]#316128 
save(adj.residuals.reffreecellmix, file="/big_data/reffree/adj.residuals_reffreecellmix_invar_filtered.Rdata")

## RUV GA
load("/big_data/reffree/adj.residuals_ruv.ga.Rdata")
adj.residuals.ruv.ga<-adj.residuals.ruv.ga[which(!(rownames(adj.residuals.ruv.ga)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.ruv.ga, file="/big_data/reffree/adj.residuals_ruv.ga_invar_filtered.Rdata")

## RUV sex
load("/big_data/reffree/adj.residuals_ruv.sex.Rdata")
adj.residuals.ruv.sex<-adj.residuals.ruv.sex[which(!(rownames(adj.residuals.ruv.sex)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.ruv.sex, file="/big_data/reffree/adj.residuals_ruv.sex_invar_filtered.Rdata")

## SVA - Supervised GA 
load("/big_data/reffree/adj.residuals_sva.sup.ga.Rdata")
adj.residuals.sva.sup.ga<-adj.residuals.sva.sup.ga[which(!(rownames(adj.residuals.sva.sup.ga)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.sva.sup.ga, file="/big_data/reffree/adj.residuals_sva.sup.ga_invar_filtered.Rdata")

## SVA - Unsupervised GA 
load("/big_data/reffree/adj.residuals_sva.unsup.ga.Rdata")
adj.residuals.sva.unsup.ga<-adj.residuals.sva.unsup.ga[which(!(rownames(adj.residuals.sva.unsup.ga)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.sva.unsup.ga, file="/big_data/reffree/adj.residuals_sva.unsup.ga_invar_filtered.Rdata")

## SVA - Supervised Sex 
load("/big_data/reffree/adj.residuals_sva.sup.sex.Rdata")
adj.residuals.sva.sup.sex<-adj.residuals.sva.sup.sex[which(!(rownames(adj.residuals.sva.sup.sex)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.sva.sup.sex, file="/big_data/reffree/adj.residuals_sva.sup.sex_invar_filtered.Rdata")

## SVA - Unsupervised Sex 
load("/big_data/reffree/adj.residuals_sva.unsup.sex.Rdata")
adj.residuals.sva.unsup.sex<-adj.residuals.sva.unsup.sex[which(!(rownames(adj.residuals.sva.unsup.sex)%in%invar_in_beta_and_independent)),] 
save(adj.residuals.sva.unsup.sex, file="/big_data/reffree/adj.residuals_sva.unsup.sex_invar_filtered.Rdata")