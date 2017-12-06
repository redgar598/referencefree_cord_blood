

###################### Figure 1 (R2 plot)
meta_cord<-read.table("~/ewas3rdround/AnalysisfileUBC.dat", sep="\t", header=T)

rownames(meta_cord)<-meta_cord$Sample_ID
facs_counts<-meta_cord[,c(7:12)]
colnames(facs_counts)<-c("NK.cells","B.cells","CD8..T.cells","CD4..T.cells","Gran","Mono")


load("~/RE_GenR/Count_like_data.Rdata")
load("~/RE_GenR/Components_GA.Rdata")
load("~/RE_GenR/Components_sex.Rdata")
load("~/RE_GenR/decon_pcs.rdata")
load("~/RE_GenR/facs_pcs.rdata")
load("~/RE_GenR/GenR_cord_deconvolution_predicted_celltypes.rdata")

load("C:/Users/040326/Documents/GenR/Count_like_data.Rdata")
load("C:/Users/040326/Documents/GenR/Components_GA.Rdata")
load("C:/Users/040326/Documents/GenR/decon_pcs.rdata")
load("C:/Users/040326/Documents/GenR/facs_pcs.rdata")
load("C:/Users/040326/Documents/GenR/GenR_cord_deconvolution_predicted_celltypes.rdata")
load("C:/Users/040326/Documents/GenR/decon_pcs.rdata")


## check all sample order
rownames(RC)<-rownames(RefFreeCounts)

rownames(ruv_GA)<-rownames(RefFreeCounts)
rownames(ruv_sex)<-rownames(RefFreeCounts)
rownames(sv_sup_sex)<-rownames(RefFreeCounts)
rownames(sv_unsup_sex)<-rownames(RefFreeCounts)
rownames(sv_sup_gestage)<-rownames(RefFreeCounts)
rownames(sv_unsup_gestage)<-rownames(RefFreeCounts)

facs_counts<-facs_counts[which(rownames(facs_counts)%in%rownames(RefFreeCounts)),]
RefFreeCounts<-RefFreeCounts[match(rownames(facs_counts), rownames(RefFreeCounts)),]
identical(rownames(RefFreeCounts), rownames(facs_counts))

RC<-RC[match(rownames(facs_counts), rownames(RC)),]
identical(rownames(RC), rownames(facs_counts))

facs_pcs<-facs_pcs[match(rownames(facs_counts), rownames(facs_pcs)),]
identical(rownames(facs_pcs), rownames(facs_counts))

decon_pcs<-decon_pcs[match(rownames(facs_counts), rownames(decon_pcs)),]
identical(rownames(decon_pcs), rownames(facs_counts))

est_cell_counts<-est_cell_counts[match(rownames(facs_counts), rownames(est_cell_counts)),]
identical(rownames(est_cell_counts), rownames(facs_counts))

ruv_GA<-ruv_GA[match(rownames(facs_counts), rownames(ruv_GA)),]
ruv_sex<-ruv_sex[match(rownames(facs_counts), rownames(ruv_sex)),]
sv_sup_sex<-sv_sup_sex[match(rownames(facs_counts), rownames(sv_sup_sex)),]
sv_unsup_sex<-sv_unsup_sex[match(rownames(facs_counts), rownames(sv_unsup_sex)),]
sv_sup_gestage<-sv_sup_gestage[match(rownames(facs_counts), rownames(sv_sup_gestage)),]
sv_unsup_gestage<-sv_unsup_gestage[match(rownames(facs_counts), rownames(sv_unsup_gestage)),]

identical(rownames(ruv_GA), rownames(facs_counts))
identical(rownames(ruv_sex), rownames(facs_counts))
identical(rownames(sv_sup_sex), rownames(facs_counts))
identical(rownames(sv_unsup_sex), rownames(facs_counts))
identical(rownames(sv_sup_gestage), rownames(facs_counts))
identical(rownames(sv_unsup_gestage), rownames(facs_counts))





### LRT for how many componenets to model each FACS cell type
library(lmtest)
## in the refactor paper they do if for many levels

## refactor on FACS

#ReFACTor	 from	 a	 linear	 model	 containing	 only	 those components	which	were	significant	(P<0.05)	under	a	nested	model	LRT

refactor_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], RC))
  colnames(countscombo)<-c("facs_count", colnames(RC))
  
  components <- colnames(RC)
  
  ## PC pval
  comp1_PC<-summary(lm(facs_count ~ 1+ PC1, data=countscombo))$coef[2,4]

  ## likelihood ratio test for refactor
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]
    }
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5)))
  
  fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
  
  formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
  formoverfit<-as.formula(paste(c("facs_count ~ 1", components[1:5]), collapse=" + "))
  
  ## now that you have the number of compenents to fit, get the r2 for the overfit and correctly fit model
  # still want the overfit model because that is realistic to refactor recommeded settings for a user with no FACS
  
  mod_refactor<-lm(formbest, data=countscombo)
  mod_refactoroverfit<-lm(formoverfit, data=countscombo)
  
  data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_refactor)$r.squared,R2_overfit=summary(mod_refactoroverfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  
})

refactor_R2<-do.call(rbind, refactor_R2)



########## RefFreeCounts

refreecellmix_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], RefFreeCounts))
  colnames(countscombo)<-c("facs_count", paste("comp",colnames(RefFreeCounts), sep=""))
  
  components <- paste("comp",colnames(RefFreeCounts), sep="")
  
  ## PC pval
  comp1_PC<-summary(lm(facs_count ~ 1+ comp1, data=countscombo))$coef[2,4]
  
  ## likelihood ratio test for refactor
  #x 1:5
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5)))
  
  fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
  
  formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
  formoverfit<-as.formula(paste(c("facs_count ~ 1", components[1:5]), collapse=" + "))
  
  ## now that you have the number of compenents to fit, get the r2 for the overfit and correctly fit model
  # still want the overfit model because that is realistic to refactor recommeded settings for a user with no FACS
  
  mod_refactor<-lm(formbest, data=countscombo)
  mod_refactoroverfit<-lm(formoverfit, data=countscombo)
  
  data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_refactor)$r.squared,R2_overfit=summary(mod_refactoroverfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  
})

refreecellmix_R2<-do.call(rbind, refreecellmix_R2)




########## FACS PCA

facsPCA_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], facs_pcs))
  colnames(countscombo)<-c("facs_count", colnames(facs_pcs))
  
  components <- colnames(facs_pcs)[1:5]
  
  ## PC pval
  comp1_PC<-summary(lm(facs_count ~ 1+ PC1, data=countscombo))$coef[2,4]
  
  ## likelihood ratio test for refactor
  #x 1:5
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5)))
  
  fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
  
  formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
  formoverfit<-as.formula(paste(c("facs_count ~ 1", components[1:5]), collapse=" + "))
  
  ## now that you have the number of compenents to fit, get the r2 for the overfit and correctly fit model
  # still want the overfit model because that is realistic to refactor recommeded settings for a user with no FACS
  
  mod_lrt<-lm(formbest, data=countscombo)
  mod_overfit<-lm(formoverfit, data=countscombo)
  
  data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_lrt)$r.squared,R2_overfit=summary(mod_overfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  
})

facs_PCA_R2<-do.call(rbind, facsPCA_R2)




########## Deconvolution counts
## this one just get the R2 fr each predicted cell type to its facs actual
est_cell_counts<-as.data.frame(est_cell_counts)
facs_count<-as.data.frame(facs_counts)

Monocytes<-summary(lm(facs_count$Mono ~ 1+ est_cell_counts$Mono))$r.squared
gran<-summary(lm(facs_count$Gran ~ 1+ est_cell_counts$Gran))$r.squared
NK.cells<-summary(lm(facs_count$NK.cells ~ 1+ est_cell_counts$NK))$r.squared
B.cells<-summary(lm(facs_count$B.cells ~ 1+ est_cell_counts$Bcell))$r.squared
CD4Tcells<-summary(lm(facs_count$CD4..T.cells ~ 1+ est_cell_counts$CD4T))$r.squared
CD8Tcells<-summary(lm(facs_count$CD8..T.cells ~ 1+ est_cell_counts$CD8T))$r.squared


deconcounts_R2<-data.frame(facs_celltype=colnames(facs_count), 
                           R2_LRT=c(Monocytes,gran,NK.cells,B.cells,CD4Tcells,CD8Tcells), 
                           R2_overfit=NA, components=NA)





########## Deconvolution PCA

deconPCA_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], decon_pcs))
  colnames(countscombo)<-c("facs_count", colnames(decon_pcs))
  
  components <- colnames(decon_pcs)[1:5]
  
  ## PC pval
  comp1_PC<-summary(lm(facs_count ~ 1+ PC1, data=countscombo))$coef[2,4]
  
  ## likelihood ratio test for refactor
  #x 1:5
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5)))
  
  fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
  
  formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
  formoverfit<-as.formula(paste(c("facs_count ~ 1", components[1:5]), collapse=" + "))
  
  ## now that you have the number of compenents to fit, get the r2 for the overfit and correctly fit model
  # still want the overfit model because that is realistic to refactor recommeded settings for a user with no FACS
  
  mod_lrt<-lm(formbest, data=countscombo)
  mod_overfit<-lm(formoverfit, data=countscombo)
  
  data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_lrt)$r.squared,R2_overfit=summary(mod_overfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  
})

deconPCA_R2<-do.call(rbind, deconPCA_R2)


########## SVA and ruv function

r2df<-function(component_matrix) {
  component_matrix<-as.data.frame(component_matrix)
  method_R2<-lapply(1:ncol(facs_counts), function(type){
    countscombo<-as.data.frame(cbind(facs_counts[,type], component_matrix))
    colnames(countscombo)<-c("facs_count", colnames(component_matrix))
    components <- colnames(component_matrix)
    
    ## PC pval
    comp1_PC<-summary(lm(facs_count ~ 1+ V1, data=countscombo))$coef[2,4]
    
    ## likelihood ratio test
    lrt<-function(x){
      form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
      form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
      mod0<-lm(form0, data=countscombo)
      mod1<-lm(form1, data=countscombo)
      lrtmod<-lrtest(mod1, mod0)
      lrtmod$`Pr(>Chisq)`[2]}
    
    lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5)))
    fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
    
    formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
    formoverfit<-as.formula(paste(c("facs_count ~ 1", components), collapse=" + "))
    mod_lrt<-lm(formbest, data=countscombo)
    mod_overfit<-lm(formoverfit, data=countscombo)
    data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_lrt)$r.squared,R2_overfit=summary(mod_overfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  })
  
  do.call(rbind, method_R2)}


ruv_GA_R2<-r2df(ruv_GA)
ruv_sex_R2<-r2df(ruv_sex)
sv_sup_sex_R2<-r2df(sv_sup_sex)
sv_unsup_sex_R2<-r2df(sv_unsup_sex)
sv_sup_gestage_R2<-r2df(sv_sup_gestage)
sv_unsup_gestage_R2<-r2df(sv_unsup_gestage)


#simplier plot
refreecellmix_R2_plt<-refreecellmix_R2[,1:2]
refreecellmix_R2_plt$method<-"RefFreeCellMix"
refactor_R2_plt<-refactor_R2[,1:2]
refactor_R2_plt$method<-"ReFACTor"
facs_PCA_R2_plt<-facs_PCA_R2[,1:2]
facs_PCA_R2_plt$method<-"FACS - PCA - Gold-Standard"
deconcounts_R2_plt<-deconcounts_R2[,1:2]
deconcounts_R2_plt$method<-"Deconvolution - Counts"
deconPCA_R2_plt<-deconPCA_R2[,1:2]
deconPCA_R2_plt$method<-"Deconvolution - PCA"
ruv_GA_R2_plt<-ruv_GA_R2[,1:2]
ruv_GA_R2_plt$method<-"RUV - GA"
ruv_sex_R2_plt<-ruv_sex_R2[,1:2]
ruv_sex_R2_plt$method<-"RUV - Sex"
sv_sup_sex_R2_plt<-sv_sup_sex_R2[,1:2]
sv_sup_sex_R2_plt$method<-"SVA - Supervised Sex"
sv_unsup_sex_R2_plt<-sv_unsup_sex_R2[,1:2]
sv_unsup_sex_R2_plt$method<-"SVA - Unsupervised Sex"
sv_sup_gestage_R2_plt<-sv_sup_gestage_R2[,1:2]
sv_sup_gestage_R2_plt$method<-"SVA - Supervised GA"
sv_unsup_gestage_R2_plt<-sv_unsup_gestage_R2[,1:2]
sv_unsup_gestage_R2_plt$method<-"SVA - Unsupervised GA"

plt_r2<-rbind(refactor_R2_plt, refreecellmix_R2_plt,facs_PCA_R2_plt,deconcounts_R2_plt,deconPCA_R2_plt,
              ruv_GA_R2_plt, ruv_sex_R2_plt,sv_sup_sex_R2_plt,sv_unsup_sex_R2_plt,sv_sup_gestage_R2_plt,sv_unsup_gestage_R2_plt)


save(plt_r2, file="~/RE_GenR/R2_plot.RData")


#################################################################################################################################
## Reverse R2
#################################################################################################################################

# refactor

R2_inverse<-function(comp_df, method_name){
  r2s<-lapply(1:ncol(comp_df), function(comp){
  data.frame(Method=method_name, comp=comp, r2=summary(lm(comp_df[,comp] ~ 1+ Mono + Gran  + NK.cells +  B.cells + CD4..T.cells + CD8..T.cells, data=facs_counts))$r.squared)})
  do.call(rbind, r2s)
}

ReFACTor<- R2_inverse(RC, "ReFACTor")
decon1<- R2_inverse(decon_pcs[,1:5], "Deconvolution - PCA")
facs2<-R2_inverse(facs_pcs[,1:5], "FACS - PCA - Gold-Standard")
reffreecellmix<-  R2_inverse(RefFreeCounts, "RefFreeCellMix")
ruvga<-R2_inverse(ruv_GA, "RUV - GA")
ruvsex<-R2_inverse(ruv_sex, "RUV - Sex")
sv1<-R2_inverse(sv_sup_gestage, "SVA - Supervised GA")
sv2<-R2_inverse(sv_sup_sex, "SVA - Supervised Sex")
sv3<-R2_inverse(sv_unsup_gestage, "SVA - Unsupervised GA")
sv4<-R2_inverse(sv_unsup_sex, "SVA - Unsupervised Sex")

## decon counts special case again
decon2<-R2_inverse(est_cell_counts, "Deconvolution - Counts")


Monocytes<-summary(lm(est_cell_counts$Mono ~ 1+ facs_count$Mono))$r.squared
gran<-summary(lm(est_cell_counts$Gran ~ 1+ facs_count$Gran))$r.squared
NK.cells<-summary(lm(est_cell_counts$NK ~ 1+ facs_count$NK.cells))$r.squared
B.cells<-summary(lm(est_cell_counts$Bcell ~ 1+ facs_count$B.cells))$r.squared
CD4Tcells<-summary(lm(est_cell_counts$CD4T ~ 1+ facs_count$CD4..T.cells))$r.squared
CD8Tcells<-summary(lm(est_cell_counts$CD8T ~ 1+ facs_count$CD8..T.cells))$r.squared


decon2<-data.frame(Method="Deconvolution - Counts", 
                   comp=1:6,
                   r2=c(Monocytes,gran,NK.cells,B.cells,CD4Tcells,CD8Tcells))

R2_inverse_all<-rbind(ReFACTor,decon1,decon2, facs2, reffreecellmix,ruvga,ruvsex, sv1,sv2,sv3,sv4)
R2_inverse_all$R2_inverse<-1-R2_inverse_all$r2

save(R2_inverse_all, file="~/RE_GenR/R2_inverse_plot.RData" )



