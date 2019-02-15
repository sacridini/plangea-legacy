# VERSION 0.2 (DRAFT 3 June 2018) CONSISTING OF 5 COMPONENTS:
# preprocessing.r, optimisation.r, run_optimisation.r, postprocessing.r, functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au

library(raster)

dir <- "inputdata_v3/"

#outdir <- "opt_results_t350_v5/"
outdir <- "opt_results_t150_v5/"

# constant to convert tonnes/ha to giga tonnes (also = peta grammes)
const_cb <- 1E-7

# constant to convert oc to $billions
const_oc <- 1E-9

# constant to convert area to millions of km2
const_area <- 1E-6

# constant to convert mean abs extinction risk delta proportion to a percent
const_bd <- 1E2


postprocess.results(dir, outdir, "scen_cb-bd", hasweights=TRUE)
postprocess.results(dir, outdir, "scen_cb-bd-oc", hasweights=TRUE)
postprocess.results(dir, outdir, "scen_cb", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_cb-oc", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_bd", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_bd-oc", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_oc", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_world", hasweights=FALSE)
for (i in 1:10) postprocess.results(dir, outdir, paste0("scen_rnd", i), hasweights=FALSE)

outdir <- "opt_results_t150_v5/"

postprocess.results(dir, outdir, "scen_cb-bd", hasweights=TRUE)
postprocess.results(dir, outdir, "scen_cb-bd-oc", hasweights=TRUE)
postprocess.results(dir, outdir, "scen_cb", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_cb-oc", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_bd", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_bd-oc", hasweights=FALSE)
postprocess.results(dir, outdir, "scen_oc", hasweights=FALSE)
for (i in 1:10) postprocess.results(dir, outdir, paste0("scen_rnd", i), hasweights=FALSE)


# COMPILE AND EXPORT ALL RESULTS 

outdir <- "opt_results_t150_v5/"

load(file=paste0(outdir, "scen_cb-bd", "_results_df.RData"))
results.df <- df
load(file=paste0(outdir, "scen_cb-bd-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_cb", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_cb-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_bd", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_bd-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
for (i in 1:10) { 
	load(file=paste0(outdir, "scen_rnd", i, "_results_df.RData"))
	results.df <- rbind(results.df, df)
}
load(file=paste0(outdir, "scen_world", "_results_df.RData"))
results.df <- rbind(results.df, df)


outdir <- "opt_results_t150_v5/"

load(file=paste0(outdir, "scen_cb-bd", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_cb-bd-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_cb", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_cb-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_bd", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_bd-oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
load(file=paste0(outdir, "scen_oc", "_results_df.RData"))
results.df <- rbind(results.df, df)
for (i in 1:10) { 
	load(file=paste0(outdir, "scen_rnd", i, "_results_df.RData"))
	results.df <- rbind(results.df, df)
}

results.df

save(results.df, file=paste0(outdir, "allscenarios_results.df.RData"))
write.csv(results.df, file=paste0(outdir, "allscenarios_results.csv"), row.names=F)

#scp 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/opt_results_t150_v4/allscen* .


### MAIN PLOT - TRAEOFF CURVES ###################################################


#load("opt_results_t150_v1/allscenarios_results.df.RData")
load(paste0(outdir, "allscenarios_results.df.RData"))

sym.cex <- 2

# colour ramp for costs
colfunc <- colorRampPalette(c("purple", "blue", "green", "yellow", "orange", "red1", "red2", "red3", "red4"))
refcols <- colfunc(100)

min.col <- 0
max.col <- 18


pdf(file="main_result_plot.pdf", width=6, height=5)
lo <- layout(matrix(c(1, 2), nrow=1), widths=c(4, 1))

par(mar=c(4, 4, 1, 1), mgp=c(2.4, 1, 0))
plot(c(0, 70), c(0, 4.5), type="n", xlab="Delta carbon sequestered (Pg CO2Eq)", ylab="Delta mean extinction risk (%)")

# random point
# points(rnd.cb, rnd.bd, col="black", bg=rnd.col, pch=21, cex=sym.cex)
# text(rnd.cb, rnd.bd, labels="V", pos=3, offset=0.7)

# IV min cost
# points(four.cb, four.bd, col="black", bg=four.col, pch=21, cex=sym.cex)
# text(four.cb, four.bd, labels="IV", pos=3, offset=0.7)

# I BAU
# points(one.cb, one.bd, col="black", bg=one.col, pch=21, cex=sym.cex)
# text(one.cb, one.bd, labels="I", pos=3, offset=0.7)

# line 1: cost effective
recs <- which(results.df$scenario == "scen_cb-bd-oc" & results.df$area.rest > 3)
lines(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], lty=2, lwd=2)
line.col <- refcols[trunc(100*((results.df$oc[recs] - min.col)/(max.col - min.col))) + 1]
points(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], col=line.col, pch=16, cex=sym.cex)
# black outline
# pid <- 1
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="III", pos=2, offset=1)
# pid <- 8
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="II", pos=2, offset=1)
# pid <- 9
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="I", pos=2, offset=1)

recs <- which(results.df$scenario == "scen_cb-bd-oc" & results.df$area.rest < 3)
lines(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], lty=2, lwd=2)
line.col <- refcols[trunc(100*((results.df$oc[recs] - min.col)/(max.col - min.col))) + 1]
points(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], col=line.col, pch=16, cex=sym.cex)


# line 2: max return
recs <- which(results.df$scenario == "scen_cb-bd" & results.df$area.rest > 3)
lines(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], lty=1, lwd=2)
line.col <- refcols[trunc(100*((results.df$oc[recs] - min.col)/(max.col - min.col))) + 1]
points(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], col=line.col, pch=16, cex=sym.cex)
# black outline
# pid <- 1
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="IV", pos=1, offset=1)
# pid <- 5
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="V", pos=4, offset=1)
# pid <- 7
# points(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], col="black", pch=1, cex=sym.cex)
# text(results.df$cb[recs[pid]], results.df$exrisk.red.mean.abs[recs[pid]], labels="VI", pos=1, offset=1)

recs <- which(results.df$scenario == "scen_cb-bd" & results.df$area.rest < 3)
lines(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], lty=1, lwd=2)
line.col <- refcols[trunc(100*((results.df$oc[recs] - min.col)/(max.col - min.col))) + 1]
points(results.df$cb[recs], results.df$exrisk.red.mean.abs[recs], col=line.col, pch=16, cex=sym.cex)


# random
s <- substr(results.df$scenario, 1, 8)
recs <- which(s == "scen_rnd" & results.df$area.rest > 3)
line.col <- refcols[trunc(100*((mean(results.df$oc[recs]) - min.col)/(max.col - min.col))) + 1]
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col=line.col, pch=16, cex=sym.cex)
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col="black", pch=1, cex=sym.cex)
text(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), labels="I", pos=1, offset=1)


# random
s <- substr(results.df$scenario, 1, 8)
recs <- which(s == "scen_rnd" & results.df$area.rest < 3)
line.col <- refcols[trunc(100*((mean(results.df$oc[recs]) - min.col)/(max.col - min.col))) + 1]
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col=line.col, pch=16, cex=sym.cex)
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col="black", pch=1, cex=sym.cex)
text(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), labels="II", pos=1, offset=1)


# oc only
s <- substr(results.df$scenario, 1, 8)
recs <- which(s == "scen_oc" & results.df$area.rest > 3)
line.col <- refcols[trunc(100*((mean(results.df$oc[recs]) - min.col)/(max.col - min.col))) + 1]
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col=line.col, pch=16, cex=sym.cex)
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col="black", pch=1, cex=sym.cex)
text(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), labels="III", pos=1, offset=1)

s <- substr(results.df$scenario, 1, 8)
recs <- which(s == "scen_oc" & results.df$area.rest < 3)
line.col <- refcols[trunc(100*((mean(results.df$oc[recs]) - min.col)/(max.col - min.col))) + 1]
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col=line.col, pch=16, cex=sym.cex)
points(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), col="black", pch=1, cex=sym.cex)
text(mean(results.df$cb[recs]), mean(results.df$exrisk.red.mean.abs[recs]), labels="IV", pos=1, offset=1)



legend("bottomright", lty=c(1, 2), legend=c("Max. Env. Frontier (Excl. Costs)", "Cost-effective"), bty="n", cex=0.7)


# plot 2 - cost legend

legend_image <- as.raster(matrix(rev(colfunc(100)), ncol=1))
par(mar=c(7, 0, 4, 3))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.5, y = seq(0,1,l=3), labels = c(2, 10, 18))
rasterImage(legend_image, 0, 0, 1,1)
mtext("Total costs (trillion USD)", side=4, line=1.2)

dev.off()













# 
# ### IGNORE EVERYTHING BELOW THIS POINT ####################################################
# 
# 
# # code for processing outputs from 
# 
# # table of return values and costs
# scen <- "scen-cb-bd"
# 
# load(file=paste0(outdir, scen, "_weightm.RData"))
# load(file=paste0(dir, "exrisk.t0.RData"))
# load(file=paste0(dir, "cb.RData"))
# load(file=paste0(dir, "ocg.RData"))
# load(file=paste0(dir, "occ.RData"))
# load(file=paste0(dir, "prop.crop.RData"))
# load(file=paste0(dir, "prop.cultg.RData"))
# oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg 
# load(file=paste0(dir, "A.RData"))
# 
# # constant to convert tonnes/ha to giga tonnes (also = peta grammes)
# const_cb <- 1E-7
# 
# # constant to convert oc to $billions
# const_oc <- 1E-9
# 
# # constant to convert area to millions of km2
# const_area <- 1E-6
# 
# # constant to convert mean abs extinction risk delta proportion to a percent
# const_bd <- 1E2
# 
# exrisk.mean <- rep(0, dim(weightm)[1])
# exrisk.sum <- rep(0, dim(weightm)[1])
# exrisk.red.mean.abs <- rep(0, dim(weightm)[1])
# exrisk.red.mean.prop <- rep(0, dim(weightm)[1])
# cb.total <- rep(0, dim(weightm)[1])
# oc.total <- rep(0, dim(weightm)[1])
# area.rest <- rep(0, dim(weightm)[1])
# 
# for (w in 1:dim(weightm)[1]){
# 
# 	load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# 	nsteps <- dim(res.exrisk)[2]
# 	exrisk.mean[w] <- mean(res.exrisk[,nsteps])
# 	exrisk.sum[w] <- sum(res.exrisk[,nsteps])
# 	exrisk.red.mean.abs[w] <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# 	#exrisk.red.mean.prop[w] <-  mean((exrisk.t0 - res.exrisk[,nsteps]) / exrisk.t0)
# 
# 	load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# 	cb.total[w] <- sum(cb * res.total.restored.pu * A * const_cb)
# 	oc.total[w] <- sum(oc * res.total.restored.pu * A * const_oc)
# 	area.rest[w] <- sum(res.total.restored.pu * A * const_area)
# 
# }
# 
# df <- data.frame(scenario=rep(scen, dim(weightm)[1]), weightcb=weightm[,1], weightbd=weightm[,2], exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# 
# # table of return values and costs
# scen <- "scen4b"
# 
# load(file=paste0(outdir, scen, "_weightm.RData"))
# load(file=paste0(dir, "exrisk.t0.RData"))
# load(file=paste0(dir, "cb.RData"))
# load(file=paste0(dir, "ocg.RData"))
# load(file=paste0(dir, "occ.RData"))
# load(file=paste0(dir, "prop.crop.RData"))
# load(file=paste0(dir, "prop.cultg.RData"))
# oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg
# load(file=paste0(dir, "A.RData"))
# 
# exrisk.mean <- rep(0, dim(weightm)[1])
# exrisk.sum <- rep(0, dim(weightm)[1])
# exrisk.red.mean.abs <- rep(0, dim(weightm)[1])
# exrisk.red.mean.prop <- rep(0, dim(weightm)[1])
# cb.total <- rep(0, dim(weightm)[1])
# oc.total <- rep(0, dim(weightm)[1])
# area.rest <- rep(0, dim(weightm)[1])
# 
# for (w in 1:dim(weightm)[1]){
# 
# 	load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# 	nsteps <- dim(res.exrisk)[2]
# 	exrisk.mean[w] <- mean(res.exrisk[,nsteps])
# 	exrisk.sum[w] <- sum(res.exrisk[,nsteps])
# 	exrisk.red.mean.abs[w] <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# 	#exrisk.red.mean.prop[w] <-  mean((exrisk.t0 - res.exrisk[,nsteps]) / exrisk.t0)
# 
# 	load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# 	cb.total[w] <- sum(cb * res.total.restored.pu * A * const_cb)
# 	oc.total[w] <- sum(oc * res.total.restored.pu * A * const_oc)
# 	area.rest[w] <- sum(res.total.restored.pu * A * const_area)
# 
# }
# 
# df <- data.frame(scenario=rep(scen, dim(weightm)[1]), weightcb=weightm[,1], weightbd=weightm[,2], exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# # load(file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
# # sum(res.total.restored.spp)
# 
# # tmp <- extinction.risk(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# # sum(tmp)
# 
# 
# 
# # CARBON ONLY
# scen <- "scen1a"
# w <- "NA"
# 
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# nsteps <- dim(res.exrisk)[2]
# exrisk.mean <- mean(res.exrisk[,nsteps])
# exrisk.sum <- sum(res.exrisk[,nsteps])
# exrisk.red.mean.abs <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# cb.total <- sum(cb * res.total.restored.pu * A * const_cb)
# oc.total <- sum(oc * res.total.restored.pu * A * const_oc)
# area.rest <- sum(res.total.restored.pu * A * const_area)
# 
# df <- data.frame(scenario=scen, weightcb=NA, weightbd=NA, exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# scen <- "scen1b"
# w <- "NA"
# 
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# nsteps <- dim(res.exrisk)[2]
# exrisk.mean <- mean(res.exrisk[,nsteps])
# exrisk.sum <- sum(res.exrisk[,nsteps])
# exrisk.red.mean.abs <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# cb.total <- sum(cb * res.total.restored.pu * A * const_cb)
# oc.total <- sum(oc * res.total.restored.pu * A * const_oc)
# area.rest <- sum(res.total.restored.pu * A * const_area)
# 
# df <- data.frame(scenario=scen, weightcb=NA, weightbd=NA, exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# 
# 
# 
# # BD ONLY
# scen <- "scen2a"
# w <- "NA"
# 
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# nsteps <- dim(res.exrisk)[2]
# exrisk.mean <- mean(res.exrisk[,nsteps])
# exrisk.sum <- sum(res.exrisk[,nsteps])
# exrisk.red.mean.abs <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# cb.total <- sum(cb * res.total.restored.pu * A * const_cb)
# oc.total <- sum(oc * res.total.restored.pu * A * const_oc)
# area.rest <- sum(res.total.restored.pu * A * const_area)
# 
# df <- data.frame(scenario=scen, weightcb=NA, weightbd=NA, exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# scen <- "scen2b"
# w <- "NA"
# 
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# nsteps <- dim(res.exrisk)[2]
# exrisk.mean <- mean(res.exrisk[,nsteps])
# exrisk.sum <- sum(res.exrisk[,nsteps])
# exrisk.red.mean.abs <-  mean(exrisk.t0 - res.exrisk[,nsteps]) * const_bd
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# cb.total <- sum(cb * res.total.restored.pu * A * const_cb)
# oc.total <- sum(oc * res.total.restored.pu * A * const_oc)
# area.rest <- sum(res.total.restored.pu * A * const_area)
# 
# df <- data.frame(scenario=scen, weightcb=NA, weightbd=NA, exrisk.mean, exrisk.sum, exrisk.red.mean.abs, cb.total, oc.total, area.rest)
# df
# save(df, file=paste0(outdir, scen, "_results_df.RData"))
# 
# 
# 
# # COMPILE AND EXPORT ALL RESULTS 
# 
# scen <- "scen1a"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- df
# 
# scen <- "scen1b"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- rbind(results.df, df)
# scen <- "scen2a"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- rbind(results.df, df)
# scen <- "scen2b"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- rbind(results.df, df)
# scen <- "scen4a"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- rbind(results.df, df)
# scen <- "scen4b"
# load(file=paste0(outdir, scen, "_results_df.RData"))
# results.df <- rbind(results.df, df)
# results.df
# 
# save(results.df, file=paste0(outdir, "allscenarios_results.df.RData"))
# write.csv(results.df, file=paste0(outdir, "allscenarios_results.csv"), row.names=F)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# scen <- "scen4b-2"
# write.csv(df, file=paste0(outdir, scen, "_results_table.csv", rownames=F))
# save(df, file=paste0(outdir, scen, "_results_table.RData"))
# 
# 
# load(file=paste0(outdir, scen, "_results_table.RData"))
# 
# 
# pdf(file="tradeoff.pdf", width=5, height=5)
# plot(df$cb.total, df$exrisk.sum, type="l", xlab="carbon", ylab="extinction risk")
# dev.off()
# 
# 
# plot(df$cb.total, type="l", xlab="carbon", ylab="extinction risk")
# 
# 
# 
# 
# # save(res.objval, file=paste0(outdir, scen, "_res.objval_w_", w, ".RData"))
# # save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# # save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# # save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# # save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
# 
# # scenario name - BD only
# scen <- "scen5-3"
# 
# load(file=paste0(dir, "spid_proc.RData"))
# load(file=paste0(dir, "spcat_proc.RData"))
# load(file=paste0(outdir, scen, "_res.exrisk_w_NA.RData"))
# 
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# dim(res.prop.restored.pu)
# 
# dim(res.exrisk)
# 
# colnames(res.exrisk) <- c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20")
# df <- data.frame(SPID=spid_proc, SPCAT=spcat_proc, exrisk.t0, res.exrisk)
# 
# write.csv(df, file="biodiversity_restoration_increments_to_35Mkm2.csv")
# 
# 
# pdf(file="bd_exrisk_profile_35Mkm2.pdf", width=5, height=5)
# plot(c(0:20)*5, 803.7516 - c(803.7516, colSums(res.exrisk)), type="l", xlab="proportion restored (of 35 M km2 max)", ylab="extinction risk reduction")
# dev.off()
# 
# 
# v <- rep(0, length(master_index))
# for (i in 1:20){
# 	recs <- which(res.prop.restored.pu[,i] > 0)
# 	v[recs] <- i
# }
# table(v)
# 
# plot.pu.map(v, fname="biodiversity_restoration_order_map4.png", bias=1)
# 
# plot.pu.map(exrisk.t0.slope, fname="biodiversity_exrisk.t0.slope_map.png", bias=1)
# 
# 
# 
# 
# 
# 
# write.pu.raster(bd, "bd_values.tif")
# 
# 
# # sunday plots
# 
# 
# load()
# 
# 
# 
# 
# dir <- "inputdata_v2/"
# 
# 
# # total terrestrial maximum potential habitat area (reference data for extinction risk calculation)
# load(file=paste0(dir, "habarea.max.RData"))
# 
# # total terrestrial habitat area at t0
# load(file=paste0(dir, "habarea.t0.RData"))
# 
# # extinction risk at t0 and slopes
# load(file=paste0(dir, "exrisk.t0.RData"))
# load(file=paste0(dir, "exrisk.t0.slope.RData"))
# 
# 
# 
# hist(habarea.max)
# 
# plot(log(habarea.t0+1), log(habarea.max+1))
# abline(coef=c(0,1))
# 
# 
# pdf(file="exrisk_distributions.pdf", width=8, height=4)
# par(mfrow=c(1,2), mar=c(4,4, 1, 1))
# hist(habarea.t0/habarea.max, main="", col="grey80", xlab="habitat area t0 / original habitat area", ylab="frequency")
# hist(exrisk.t0, main="", col="grey80", xlab="extinction risk", ylab="frequency")
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# 
# # the species indices with reference to the master_index
# load(file=paste0(dir, "species_index_list_proc.RData"))
# 
# # load(file=paste0(dir, "upperbound.RData"))
# 
# # proportion of habitat in each planning unit cell
# # load(file=paste0(dir, "prop.hab.pu.RData"))
# 
# # species habitat associations
# load(file=paste0(dir, "sphabm_proc.RData"))
# 
# 
# # load other covariate data
# load(file=paste0(dir, "cb.RData"))
# load(file=paste0(dir, "ocg.RData"))
# load(file=paste0(dir, "occ.RData"))
# load(file=paste0(dir, "cntry.RData"))
# load(file=paste0(dir, "A.RData"))
# 
# # proportion of ag and pasture
# load(file=paste0(dir, "prop.crop.RData"))
# load(file=paste0(dir, "prop.cultg.RData"))
# 
# # proportions of 5 nat veg types when restoration occurs
# load(file=paste0(dir, "prop.restore.RData"))
# 
# # the proportion of nat veg habitat in each PU at start
# load(file=paste0(dir, "prop.hab.pu.RData"))
# 
# # data for speeding up BD calculations
# load(file=paste0(dir, "usphab_proc.RData"))
# load(file=paste0(dir, "usphab_index.RData"))
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # maps
# load(file=paste0(dir, "master_index.RData"))
# load(file=paste0(dir, "terrestrial_index.RData"))
# 
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=1)
# cols <- c("grey90", colfunc(50))
# 
# 
# load(file="terrestrial_lands.RData")
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# 
# r.restprop <- r.terr
# values(r.restprop) <- NA
# values(r.restprop)[terrestrial_index] <- 0
# values(r.restprop)[master_index] <- res.total.restored.pu
# 
# png(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".png"), width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.restprop, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# load(r.amph, file="r.amph.RData")
# load(r.mamm, file="r.mamm.RData")
# 
# 
# 
# 
# # plot input base maps
# 
# # for now, calculate oc as the weighted average of occ and ocg:
# 
# 
# # load other covariate data
# load(file=paste0(dir, "cb.RData"))
# load(file=paste0(dir, "ocg.RData"))
# load(file=paste0(dir, "occ.RData"))
# 
# # proportion of ag and pasture
# load(file=paste0(dir, "prop.crop.RData"))
# load(file=paste0(dir, "prop.cultg.RData"))
# 
# 
# oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg 
# 
# 
# 
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=1)
# cols <- c("grey90", colfunc(100))
# cols <- colfunc(100)
# 
# # carbon
# r <- r.terr
# values(r) <- NA
# values(r)[terrestrial_index] <- 0
# values(r)[master_index] <- cb
# png(file=paste0(dir, "base_carbon.png"), width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# # carbon/oppcost
# r <- r.terr
# values(r) <- NA
# # values(r)[terrestrial_index] <- 0
# values(r)[master_index] <- cb / oc
# png(file=paste0(dir, "base_carbon_div_oc.png"), width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("white", "grey90"), axes=F, box=F)
# plot(r, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# 
# 
# 
# # map of all terrestrial cells, and subset of cells available for restoration
# 
# 
# 
# 
# # create raster of all terrestrial cells
# 
# # r.terr <- r.er
# # recs <- which(values(r.terr) > 0)
# # values(r.terr)[recs] <- rep(1, length(recs))
# # save(r.terr, file="terrestrial_lands.RData")
# load(file="terrestrial_lands.RData")
# # writeRaster(r.terr, filename="./terrestrial_lands.tif", dataType="INT1S")
# 
# plot(r.terr, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey50"))
# 
# 
# # cells we evaluate
# r.avail <- r.terr
# 
# load(file=paste0(dir, "master_index.RData"))
# r.avail[master_index] <- 2
# 
# save(r.avail, file="restorable_lands.RData")
# load(file="restorable_lands.RData")
# 
# 
# png(file="restorable_lands.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.avail, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey50", "green"), axes=F, box=F, legend=F)
# dev.off()
# 
# 
# r.restprop <- r.terr
# values(r.restprop) <- 0
# values(r.restprop)[master_index] <- upperbound
# png(file="restoration_proportion.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# #plot(r.terr, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col="grey90", axes=F, box=F)
# plot(r.restprop, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# 
# 
# r.restprop <- r.terr
# values(r.restprop) <- 0
# values(r.restprop)[master_index] <- upperbound
# png(file="restoration_proportion.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# #plot(r.terr, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col="grey90", axes=F, box=F)
# plot(r.restprop, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# 
# 
# 
# 
# 
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=3)
# cols <- colfunc(50)
# 
# png(file="map_base_opp_cost_grass.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.ocg, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=3)
# cols <- colfunc(50)
# 
# png(file="map_base_opp_cost_crop.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.occ, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=3)
# cols <- colfunc(50)
# 
# png(file="map_base_opp_cost_avg.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.oc, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F)
# dev.off()
# 
# 
# 
# 
# 
# # plots of various z values
# 
# zs <- seq(0.1, 0.4, 0.05)
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=1)
# cols <- colfunc(length(zs))
# x <- seq(0, 10000, 10)
# 
# pdf(file="extinction_risk_z_values.pdf", width=5, height=5)
# par(mar=c(4.4, 4.4, 1, 1), mgp=c(3.2, 1, 0))
# plot(x, extinction.risk(x, Amax=10000, z=0.25), type="n", xlab="habitat area", ylab="extinction risk")
# for (i in 1:length(zs)){
# 	lines(x, extinction.risk(x, Amax=10000, z=zs[i]), col=cols[i])
# }
# legend("topright", lwd=1, col=cols, legend=zs)
# dev.off()
# 
# 
# 
# # renato's model
# 
# dados<-read.table("Plants_Forest_areas.txt", header=T, sep="\t")
# head(dados)
# 
# y <- dados[,12]
# summary(y)
# summary(dados[,9])
# x <- dados[,9]
# 
# plot(x, y)
# 
# y2 <- log(y)
# plot(x, y2)
# 
# m1 <- glm(y2 ~ x)
# summary(m1)
# 
# px <- seq(0, 80, 0.2)
# py <- predict(m1, newdata=data.frame(x=px), se.fit=T)
# 
# plot(x, y)
# lines(px, exp(py), col="red", lwd=2)
# 
# 
# 
# plot(log(x), y)
# 
# y <- y + 0.5  # <-------------- why?
# 
# m1 <- glm(y~x)
# summary(m1)
# abline(coef=coef(m1))
# px <- seq(0, 80, 10)
# py <- predict(m1, newdata=data.frame(x=px), se.fit=T)
# lines(px, py$fit + 1.96*py$se.fit, col="blue", lty=2)
# lines(px, py$fit - 1.96*py$se.fit, col="blue", lty=2)
# 
# x2 <- log(x+1)
# m2 <- glm(y ~ x2, Gamma("identity"))
# summary(m2)
# 
# 
# 
# res1 <- res.prop.restored.pu
# res1a <- rowSums(res1)
# 
# scen <- "scen_cb-bd"
# w <- 1
# 
# load(file=paste0(outdir, scen, "_res.objval_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# 
# 
# scen <- "scen_cb"
# w <- "NA"
# load(file=paste0(outdir, scen, "_delta.hab.pu_w_", w, ".RData"))
# d1 <- delta.hab.pu
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# sum(res.exrisk[,1])
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# pr1 <- res.prop.restored.pu
# load(file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# ar1 <- res.area.restored.spp
# 
# scen <- "scen_cb-bd"
# w <- 1
# load(file=paste0(outdir, scen, "_delta.hab.pu_w_", w, ".RData"))
# d2 <- delta.hab.pu
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# sum(res.exrisk[,5])
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# pr2 <- res.prop.restored.pu
# pr2s <- rowSums(pr2)
# load(file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# ar2 <- res.area.restored.spp
# 
# sum(res.total.restored.pu)
# colSums(res.exrisk)
# colSums(res.area.restored.spp)
# colSums(res.prop.restored.pu)
# 
# scen <- "scen_cb"
# w <- "NA"
# load(file=paste0(outdir, scen, "_res.objval_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# res1 <- res.prop.restored.pu
# 
# 
# scen <- "scen_cb-bd"
# w <- 1
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# 
# delta.hab.pu <- rowSums(res.prop.restored.pu)
# delta.vegtype.pu <- prop.restore * delta.hab.pu
# delta.hab.spp <- rep(0, ns)
# for (i in 1:ns){
# 	for(j in 1:5){
# 		if (sphabm_proc[i,j] == 1){
# 			delta.hab.spp[i] <- delta.hab.spp[i] + sum(delta.vegtype.pu[species_index_list_proc[[i]], j], na.rm=T)
# 		}
# 	}
# }
# delta.hab.spp <- delta.hab.spp * A
# res.total.restored.spp <- delta.hab.spp
# exrisk.slope <- calc.extinction.slope(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# res.exrisk <- extinction.risk(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# sum(res.exrisk)
# res1 <- delta.hab.pu
# er1 <- res.exrisk[,5]
# 
# scen <- "scen_cb"
# w <- "NA"
# load(file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# load(file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
# 
# delta.hab.pu <- res.prop.restored.pu[,1]
# delta.vegtype.pu <- prop.restore * delta.hab.pu
# delta.hab.spp <- rep(0, ns)
# for (i in 1:ns){
# 	for(j in 1:5){
# 		if (sphabm_proc[i,j] == 1){
# 			delta.hab.spp[i] <- delta.hab.spp[i] + sum(delta.vegtype.pu[species_index_list_proc[[i]], j], na.rm=T)
# 		}
# 	}
# }
# delta.hab.spp <- delta.hab.spp * A
# res.total.restored.spp <- delta.hab.spp
# exrisk.slope <- calc.extinction.slope(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# res.exrisk <- extinction.risk(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# sum(res.exrisk)
# res2 <- delta.hab.pu
# er2 <- res.exrisk[,1]
# 
# recs <- which(!res1 == res2)
# length(recs)
# sum(res1-res2)
# 
# 
# 
# 
# # intentionally mis-apply the calc.bd function to map extinction risk at time 0 for bernardo
# 
# tmp <- calc.bd(exrisk.t0)
# # plot.pu.map(tmp, fname="extinctionrisk_t0_map.png")
# colfunc <- colorRampPalette(c("blue", "green", "yellow", "orange", "red"), bias=1)
# cols <- colfunc(100)
# r <- r.terr
# values(r) <- NA
# values(r)[master_index] <- tmp
# png(file="extinctionrisk_t0_map.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("white", "grey90"), axes=F, box=F, legend=F)
# plot(r, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# 
# r <- r.er
# values(r) <- rep(NA, ncell(r))
# values(r)[master_index] <- tmp
# writeRaster(r, file="exrisk_t0.tif")
# 
# 
# 
