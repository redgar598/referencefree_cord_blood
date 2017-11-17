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


##########################################################################################################################

## load in testing data
library(methylumi)
load("/big_data/redgar/cordblood/2017.10.17_WB_BMIQnorm_combat.rdata")
beta<-betas(WB.bat)

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
invar_in_beta_and_independent<-intersect(y$CpG, rownames(Invariable_in_beta)) #109841/114204 (96.1%)
Betas_variable<-beta[which(!(rownames(beta)%in%invar_in_beta_and_independent)),]#325266 
beta<-Betas_variable 

## need mval for reference free
Mval<-function(beta) log2(beta/(1-beta))
mval = apply(as.data.frame(beta), 1, Mval) # need mvalues for combat
mval = as.data.frame(mval)
mval = t(mval)

## impute missing mval with row median
NA2med <- function(x) replace(x, is.na(x), median(x, na.rm = TRUE))
mval_complete <- t(apply(mval,1, NA2med))
beta_complete <- t(apply(beta,1, NA2med))

## load meta data
meta_cord<-read.csv("/big_data/redgar/cordblood/24_wholebloods_samplesheet.csv")
meta_cord<-meta_cord[match(colnames(beta), meta_cord$X),]
identical(colnames(beta), as.character(meta_cord$X))

              
################# RefFreeCellMix
library(RefFreeEWAS)
#Kˆ=3K^=3 for the cord blood data set BL-as, larger for the other three peripheral blood datasets
# We performed reference-free deconvolution with the method RefFreeCellMix by Houseman et al. [42] using the R package RefFreeEWAS. In accordance with the original publication of the method [42], we applied it to the 20,000 most variable CpG positions from the methylation matrix, unless the total number of rows was less, in which case we used the full matrix. In the former case, we used the available option to obtain the estimates of the methylation components for all CpGs as the final step of the deconvolution procedure (supplying the complete data matrix as argument Yfinal).

testArray1 <- RefFreeCellMix(as.matrix(beta_complete),K=5,iters=5,Yfinal=NULL)#verbose=T
RefFreeCounts<-as.data.frame(testArray1$Omega)


################### ReFACTor
#https://github.com/cozygene/refactor/tree/master/R
source("/big_data/redgar/cordblood/refactor.R")
        beta_df<-as.data.frame(beta_complete)
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
svobj = sva(as.matrix(mval_complete),mod,mod0,n.sv=5)
sv_unsup_gestage<-svobj$sv

################### SVA supervised GA
#Surrogate variable analysis: https://www.bioconductor.org/packages/release/bioc/html/sva.html
library(sva)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs.Rdata")
mval_diff_celltype<-mval_complete[which(rownames(mval_complete)%in%diff_cpgs$CpG),]

#Null model matrix must be nested in the full model matrix
mod = model.matrix(~meta_cord$GA)
mod0 = model.matrix(~1, data.frame(meta_cord$GA))

## surrogates
svobj = sva(as.matrix(mval_diff_celltype),mod,mod0,n.sv=5)
sv_sup_gestage<-svobj$sv

################### RUV
library(missMethyl)
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs.Rdata")

# which probes define cell type
ctl<-(rownames(mval_complete)%in%diff_cpgs$CpG)

fit = RUVfit(data=mval_complete, design=meta_cord$GA,  ctl=ctl, k=5, method=c("ruv4"))
ruv_GA<-t(fit$W)

save(sv_unsup_gestage, sv_sup_gestage, ruv_GA, file="/big_data/redgar/cordblood/Components_GA.Rdata")




##############################################################################################################################
##############################################################################################################################
load("/big_data/redgar/cordblood/Count_like_data.Rdata")
################## Corrected Betas reffreecellmix
colnames(RefFreeCounts)<-paste("comp", colnames(RefFreeCounts), sep="")

betas.lm <- apply(beta_complete, 1, function(x){
  components <- RefFreeCounts[colnames(beta_complete),]
  lm(x~comp1+comp2+comp3+comp4+comp5+0,data=components) 
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



##################### Corrected Betas refactor
rownames(RC)<-colnames(beta_complete)

betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~PC1+PC2+PC3+PC4+PC5+0,data=RC) 
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


##################### Corrected Betas SVA unsup GA
sv_unsup_ga<-as.data.frame(sv_unsup_gestage)
rownames(sv_unsup_ga)<-colnames(beta_complete)
betas.lm <- apply(beta_complete, 1, function(x){
  lm(x~V1+V2+V3+V4+V5+0,data=sv_unsup_ga) 
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
  lm(x~V1+V2+V3+V4+V5+0,data=sv_sup_ga) 
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
  lm(x~V1+V2+V3+V4+V5+0,data=ruv_GA) 
})
residuals <- t(sapply(betas.lm, function(x)residuals(summary(x))))
colnames(residuals) <- colnames(beta_complete) # re-name residuals columns with sample names
adj.residuals.ruv.ga <- residuals+matrix(apply(beta_complete, 1, mean), nrow = nrow(residuals), ncol = ncol(residuals))
save(adj.residuals.ruv.ga,file="/big_data/redgar/cordblood/adj.residuals_ruv.ga.Rdata")

rm(adj.residuals.ruv.ga)
rm(betas.lm)
gc()