load("/big_data/reffree/Cord_blood_celltype_ttest_statistics.Rdata")

head(probes_tstats[[1]][[1]])


diff_cpgs<-lapply(1:length(probes_tstats), function(x){
  df<-probes_tstats[[x]][[1]]
  top<-df[which(df$p.value<1e-08),]
  top<-top[rev(order(abs(top$dm))),]
  data.frame(CpG=rownames(top[1:100,]), celltype=names(probes_tstats)[x])
  })

diff_cpgs<-do.call(rbind,diff_cpgs)

save(diff_cpgs, file="/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs.Rdata")
load("/big_data/redgar/cordblood/Cord_blood_celltype_Diff_CpGs.Rdata")

## Meg sent list of 700 probes
load("/big_data/reffree/cord_cell_type_specific_probes.rdata")

length(intersect(as.character(probes), diff_cpgs$CpG))

## yay the probes are the same