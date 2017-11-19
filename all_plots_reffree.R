#define consistent color scheme for methods
library(ggsci)
myColors <- c("gold","goldenrod",
              "#2171b5","#6baed6","#6baed6",
              "#9467BDFF","#D62728FF",
              "#238443","#78c679","#addd8e","#d9f0a3",
              "#dd3497","#fa9fb5","grey")
show_col(myColors)

color_possibilities<-c("FACS - PCA - Gold-Standard","FACS - Drop One Cell Type",
                       "Deconvolution - PCA","Deconvolution - Drop One Cell Type","Deconvolution - Counts",
                       "ReFACTor","RefFreeCellMix",
                       "SVA - Supervised GA","SVA - Unsupervised GA",
                       "SVA - Supervised Sex","SVA - Unsupervised Sex",
                       "RUV - GA", "RUV - Sex",
                       "Uncorrected")

names(myColors) <- color_possibilities
fillscale <- scale_fill_manual(name="Method",
                               values = myColors, drop = T)







###################### Figure 1 (R2 plot)
facs_counts<-read.csv("/big_data/redgar/cordblood/Updated_counts_WB_validation.csv")

load("/big_data/redgar/cordblood/Count_like_data.Rdata")
load("/big_data/redgar/cordblood/Components_sex.Rdata")
load("/big_data/redgar/cordblood/Components_GA.Rdata")
load("/big_data/reffree/Decon_pcs.rdata")
load("/big_data/reffree/FACS_pcs.rdata")
load("/big_data/reffree/Louie_predicted_WB_celltypes_Oct27_use_Louies_sorted.rdata")

## check all sample order
rownames(facs_counts)<-facs_counts$X
facs_counts$X<-NULL

rownames(RC)<-rownames(RefFreeCounts)

rownames(ruv_GA)<-rownames(RefFreeCounts)
rownames(ruv_sex)<-rownames(RefFreeCounts)
rownames(sv_sup_sex)<-rownames(RefFreeCounts)
rownames(sv_unsup_sex)<-rownames(RefFreeCounts)
rownames(sv_sup_gestage)<-rownames(RefFreeCounts)
rownames(sv_unsup_gestage)<-rownames(RefFreeCounts)



RefFreeCounts<-RefFreeCounts[match(rownames(facs_counts), rownames(RefFreeCounts)),]
identical(rownames(RefFreeCounts), rownames(facs_counts))

RC<-RC[match(rownames(facs_counts), rownames(RC)),]
identical(rownames(RC), rownames(facs_counts))

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
est_cell_counts<-as.data.frame(est_cell_counts)
facs_count<-as.data.frame(facs_counts)
deconcounts_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], est_cell_counts))
  colnames(countscombo)<-c("facs_count", colnames(est_cell_counts))
  
  components <- colnames(est_cell_counts)
  
  ## PC pval
  comp1_PC<-summary(lm(facs_count ~ 1+ Gran, data=countscombo))$coef[2,4]
  
  ## likelihood ratio test for refactor
  #x 1:5
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(comp1_PC,lrt(2),lrt(3),lrt(4),lrt(5),lrt(6),lrt(7)))
  
  fitcomponenets<-as.character(lrtp_comp$component[which(lrtp_comp$lrtpval<0.05)])
  
  formbest<-as.formula(paste(c("facs_count ~ 1", fitcomponenets), collapse=" + "))
  formoverfit<-as.formula(paste(c("facs_count ~ 1", components), collapse=" + "))
  
  ## now that you have the number of compenents to fit, get the r2 for the overfit and correctly fit model
  # still want the overfit model because that is realistic to refactor recommeded settings for a user with no FACS
  
  mod_lrt<-lm(formbest, data=countscombo)
  mod_overfit<-lm(formoverfit, data=countscombo)
  
  data.frame(facs_celltype=colnames(facs_counts)[type], R2_LRT=summary(mod_lrt)$r.squared,R2_overfit=summary(mod_overfit)$r.squared, components=paste(fitcomponenets, collapse=", "))
  
})

deconcounts_R2<-do.call(rbind, deconcounts_R2)







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



# ####### plot it
# refreecellmix_R2_plt<-melt(refreecellmix_R2[,1:3])
# refreecellmix_R2_plt$method<-"reffreecellmix"
# refactor_R2_plt<-melt(refactor_R2[,1:3])
# refactor_R2_plt$method<-"refactor"
# 
# plt_r2<-rbind(refactor_R2_plt, refreecellmix_R2_plt)
# plt_r2$Method_R2<-paste(plt_r2$method, plt_r2$variable)
# levels(plt_r2$variable)<-c("Model with only LRT Components","Model with all 6 Components")
# 
# ggplot(plt_r2, aes(facs_celltype, value, color=method))+geom_point()+theme_bw()+
#   scale_color_manual(values=c("#b2182b","#2166ac"))+facet_wrap(~variable)+
#   ylab("R2 (FACS Cell Type Variance Accounted for)")+xlab("FACS Cell Type")

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

ggplot(plt_r2, aes(facs_celltype, R2_LRT, fill=method))+geom_point(shape=21, color="black", size=2, position=position_jitter(width=0.2))+theme_bw()+
  fillscale+
  ylab("R2 (FACS Cell Type Variance Accounted for)")+xlab("FACS Cell Type")+
  theme(text=element_text(size=12),
        axis.title=element_text(size=14))
  

ggsave("/big_data/redgar/cordblood/figures/Figure1R2.pdf", width = 8, height = 7, units = "in")
