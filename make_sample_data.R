load("/big_data/reffree/WB_betas_BMIQ_combat_together.rdata")
beta<-as.data.frame(validation_betas.combat)


load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs_Nov21.Rdata")

sample_beta<-beta[c(which(rownames(beta)%in%diff_cpgs$CpG), sample(1:nrow(beta), 1000)), ]
save(sample_beta, file="sample_beta_forCCE.RData")
mval_diff_celltype<-mval_variable[which(rownames(mval_variable)%in%diff_cpgs$CpG),]

## load meta data
meta_cord<-read.csv("/big_data/redgar/cordblood/24_wholebloods_samplesheet.csv")
meta_cord<-meta_cord[match(colnames(beta), meta_cord$X),]
identical(colnames(beta_variable), as.character(meta_cord$X))
