load("~/ewas3rdround/adj.residuals.sva.unsup.sex_predicted_cell_types.Rdata")
sva.unsup.sex<- est_cell_counts
load("~/ewas3rdround/adj.residuals.sva.unsup.ga_predicted_cell_types.Rdata")
sva.unsup.ga<- est_cell_counts
load("~/ewas3rdround/adj.residuals.sva.sup.sex_predicted_cell_types.Rdata")
sva.sup.sex<- est_cell_counts
load("~/ewas3rdround/adj.residuals.sva.sup.ga_predicted_cell_types.Rdata")
sva.sup.ga<- est_cell_counts
load("~/ewas3rdround/adj.residuals.ruv.sex_predicted_cell_types.Rdata")
ruv.sex<- est_cell_counts
load("~/ewas3rdround/adj.residuals.ruv.ga_predicted_cell_types.Rdata")
ruv.ga<- est_cell_counts
load("~/ewas3rdround/adj.residuals.reffreecellmix_predicted_cell_types.Rdata")
cellmix<- est_cell_counts
load("~/ewas3rdround/adj.residuals.refactor_predicted_cell_types.Rdata")
refactor<- est_cell_counts
load("~/ewas3rdround/facs_predicted_cell_types.Rdata")
facs.counts<- est_cell_counts
load("~/ewas3rdround/facs_pcs_predicted_cell_types.Rdata")
facs.pcs.counts<- est_cell_counts
load("~/ewas3rdround/decon_predicted_cell_types.Rdata")
decon.counts<- est_cell_counts
load("~/ewas3rdround/decon_pcs_predicted_cell_types.Rdata")
decon.pcs.counts<- est_cell_counts

load("~/ewas3rdround/GenR_cord_deconvolution_predicted_celltypes.rdata")
uncorrected<- est_cell_counts


mads<- data.frame(FACS_counts=apply(facs.counts, 2, mad),
                  FACS_pcs=apply(facs.pcs.counts, 2, mad),
                  Decon_counts=apply(decon.counts, 2, mad),
                  Decon_pcs=apply(decon.pcs.counts, 2, mad),
                  sva.unsup.sex= apply(sva.unsup.sex, 2, mad),
                  sva.unsup.ga= apply(sva.unsup.ga, 2, mad),
                  sva.sup.sex= apply(sva.sup.sex, 2, mad),
                  sva.sup.ga= apply(sva.sup.ga, 2, mad),
                  ruv.sex= apply(ruv.sex, 2, mad),
                  ruv.ga= apply(ruv.ga, 2, mad),
                  cellmix= apply(cellmix, 2, mad),
                  refactor=apply(refactor, 2, mad),
                  uncorrected=apply(uncorrected, 2, mad)
                  
)

mads$celltypes<- rownames(mads)
mads.melt<- melt(mads, id="celltypes")
save(mads.melt, file="~/ewas3rdround/mads_melted.RData")
