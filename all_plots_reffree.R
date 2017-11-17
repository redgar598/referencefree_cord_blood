#define consistent color scheme for methods
library(ggsci)
myColors <- c("gold","goldenrod",
              "#2171b5","#6baed6",
              "#9467BDFF","#D62728FF",
              "#238443","#78c679","#addd8e","#d9f0a3",
              "#dd3497","#fa9fb5")
show_col(myColors)

color_possibilities<-c("FACS - PCA - Gold-Standard","FACS - Drop One Cell Type",
                       "Deconvolution - PCA","Deconvolution - Drop One Cell Type",
                       "ReFACTor","RefFreeCellMix",
                       "SVA - Supervised GA","SVA - Unsupervised GA",
                       "SVA - Supervised Sex","SVA - Unsupervised Sex",
                       "RUV - GA", "RUV - Sex")

names(myColors) <- color_possibilities
fillscale <- scale_fill_manual(name="Method",
                               values = myColors, drop = T)


### Figure 1 (R2 plot)
load("/big_data/redgar/cordblood/Count_like_data.Rdata")
facs_counts<-read.csv("/big_data/redgar/cordblood/Updated_counts_WB_validation.csv")

rownames(RC)<-rownames(RefFreeCounts)
facs_counts<-facs_counts[match(rownames(RC), facs_counts$X),]

identical(rownames(RC), as.character(facs_counts$X))

rownames(facs_counts)<-facs_counts$X
facs_counts$X<-NULL

colnames(RefFreeCounts)<-paste("comp", colnames(RefFreeCounts), sep="")


### LRT for how many componenets to model each FACS cell type
library(lmtest)
## in the refactor paper they do if for many levels

## refactor on FACS

#ReFACTor	 from	 a	 linear	 model	 containing	 only	 those components	which	were	significant	(P<0.05)	under	a	nested	model	LRT

refactor_R2<-lapply(1:ncol(facs_counts), function(type){
  countscombo<-as.data.frame(cbind(facs_counts[,type], RC))
  colnames(countscombo)<-c("facs_count", colnames(RC))
  
  components <- colnames(RC)
  
  ## likelihood ratio test for refactor
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(lrt(1),lrt(2),lrt(3),lrt(4),lrt(5)))
  
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
  
  ## likelihood ratio test for refactor
  #x 1:6
  lrt<-function(x){
    form0<-as.formula(paste(c("facs_count ~ 1", components[0:(x-1)]), collapse=" + "))
    form1<-as.formula(paste(c("facs_count ~ 1", components[0:x]), collapse=" + "))
    
    mod0<-lm(form0, data=countscombo)
    mod1<-lm(form1, data=countscombo)
    
    lrtmod<-lrtest(mod1, mod0)
    lrtmod$`Pr(>Chisq)`[2]}
  
  
  lrtp_comp<-data.frame(component=components, lrtpval=c(lrt(1),lrt(2),lrt(3),lrt(4),lrt(5)))
  
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

plt_r2<-rbind(refactor_R2_plt, refreecellmix_R2_plt)

ggplot(plt_r2, aes(facs_celltype, R2_LRT, fill=method))+geom_point(shape=21, color="black", size=3)+theme_bw()+
  fillscale+
  ylab("R2 (FACS Cell Type Variance Accounted for)")+xlab("FACS Cell Type")+
  theme(text=element_text(size=12),
        axis.title=element_text(size=14))
  

ggsave("/big_data/redgar/cordblood/figures/Figure1R2.pdf", width = 8, height = 7, units = "in")
