# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

# This code is memory intensive, so best to start with a fresh instance of R,
# and load as few libraries and objects as needed

# This code was modified from its previous version in order to redirect the
# processing flow at the end of it to wrap_scenarios.r, instead of
# run_scenarios.r - that is because wrap_scenarios.r defines the scenario
# parameters and the range of those to be included. Those ranges are determined
# by each scenario-block run, defined in wrap_scenario.r (from which this
# script is called).

# The code in here was also modified to remove the definition of target areas
# to restore from it. The target areas are now controlled from
# wrap_optimisation.r (one level above in the processing flow), with values
# defined in wrap_scenario.r (one level further in the processing flow).

#scp functions.r 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
#scp run_scenarios.r 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
#scp run_optimisation.r 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/

#rm(list=ls())

#library(gurobi)
library(Rsymphony)
library(Matrix)

home.dir = '~/Documents/IIS_PROJECTS/plangea-legacy/'
setwd(home.dir)

source("functions.r")

# set the input directory from which to load all the data for the optimisation
dir <- 'inputdata_v8/'

#scp ~/Documents/proj/iis_global_spp/inputdata_v1/species_index_list.RData 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/inputdata_v1/

# the species indices with reference to the master_index
load(file=paste0(dir, "species_index_list_proc.RData"))

# load(file=paste0(dir, "upperbound.RData"))

# total terrestrial maximum potential habitat area (reference data for extinction risk calculation)
load(file=paste0(dir, "habarea.max.RData"))

# total terrestrial habitat area at t0
load(file=paste0(dir, "habarea.t0.RData"))

# proportion of habitat in each planning unit cell
# load(file=paste0(dir, "prop.hab.pu.RData"))

# species habitat associations
load(file=paste0(dir, "sphabm_proc.RData"))


# extinction risk at t0 and slopes
load(file=paste0(dir, "exrisk.t0.RData"))
load(file=paste0(dir, "exrisk.t0.slope.RData"))

# load other covariate data
load(file=paste0(dir, "cb.RData"))

# must prevent negative values in the carbon:
cb <- cb + abs(min(cb, na.rm=T)) + 1
summary(cb)

# Loading Opportunity cost values
load(file=paste0(dir, "ocg.RData"))
load(file=paste0(dir, "occ.RData"))
load(file=paste0(dir, "occ.diff-DR.RData"))
load(file=paste0(dir, "ocg.diff-DR.RData"))
load(file=paste0(dir, "occ.diff-Pr.RData"))
load(file=paste0(dir, "ocg.diff-Pr.RData"))

load(file=paste0(dir, "cntry.RData"))
load(file=paste0(dir, "A.RData"))

# Load standard deviation of carbon
load(file=paste0(dir, "cb-sd.RData"))

# proportion of ag and pasture
load(file=paste0(dir, "prop.crop.RData"))
load(file=paste0(dir, "prop.cultg.RData"))

# proportions of 5 nat veg types when restoration occurs
load(file=paste0(dir, "prop.restore.RData"))

# the proportion of nat veg habitat in each PU at start
load(file=paste0(dir, "prop.hab.pu.RData"))

# data for speeding up BD calculations
load(file=paste0(dir, "usphab_proc.RData"))
load(file=paste0(dir, "usphab_index.RData"))

# Species basic data.frame
load(file=paste0(dir, 'spid_proc.RData'))
load(file=paste0(dir, 'spcat_proc.RData'))
load(file=paste0(dir, 'spp_taxon_proc.RData'))
load(file=paste0(dir, 'habarea.t0.RData'))
load(file=paste0(dir, 'habarea.max.RData'))
load(file=paste0(dir, 'exrisk.t0.RData'))

base.spp.list = data.frame('Species ID'=spid_proc,
                           'Category'=spcat_proc, 'Taxon'=spp_taxon_proc)

# Master index
load(file=paste0(dir, "master_index.RData"))

# total area available for restoration (100%) expressed as proportions of cell areas
total.restoration.area <- sum(prop.crop + prop.cultg) * A
# millions of km to restore:
total.restoration.area/1E6

# set the number of planning units and number of species
ns <- length(habarea.max)
np <- length(cb)


###  SCALING OF VARIABLES
# The relative weights between CB and BD are heavily influenced by the magnitude of those variables as well as opportunity cost, and it helps to ensure they are on roughly the same scale, so we set a global scalar for each of these and it should not change among scenarios
# Hence the weights may need adjustment among scenarios. The sole purpose of the weights is just to effectively describe the tradeoff curve, so it is irrelevant if the weights vary among scenarios

# carbon:
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# -2036.31    24.13    48.91    58.39    80.97   575.65 

g_scalar_cb <- 1E-3

# biodiversity:
# this will need to be adjusted when birds are added:
g_scalar_bd <- 1E2
#g_scalar_bd <- 1E-1  # <----- this was for testing purposes

# opportunity cost:
# > summary(ocg)
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#    193.6    371.3    530.0    962.6    870.1 121429.8
# > summary(occ)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1928    3771    5299    6075    6635  100516
g_scalar_oc <- 1E-3

# we need to be very careful about numerical issues arising from large values
# like target areas in the constraits too
g_scalar_area <- 1E-4

g_diagnostic_plots <- TRUE
g_diagnostic_plots <- FALSE

# flag to help with monitoring and error detection
g_fatal_error <- FALSE

# the cb and oc rescaling can occur immediately:
summary(cb)
cb <- cb * g_scalar_cb
summary(cb)

cb.sd = cb.sd * g_scalar_cb

# for now, calculate oc as the weighted average of occ and ocg:
oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg 
summary(oc)

# Rescaling OC to be input to solver
oc <- oc * g_scalar_oc
summary(oc)


# OC standard deviation

# SD on opportunity costs propagated from uncertainty on the Discount Rate

occ.sd = sqrt(((occ.diff.Pr^2) * ((PrC.relative.var*median(occ))^2)) + ((occ.diff.DR^2) * (DR.sd^2)))
summary(occ.sd)

ocg.sd = sqrt(((ocg.diff.Pr^2) * ((PrG.relative.var*median(ocg))^2)) + ((ocg.diff.DR^2) * (DR.sd^2)))
summary(ocg.sd)

oc.sd = sqrt(
  ((prop.crop / (prop.crop + prop.cultg))^2) * (occ.sd^2) +
    ((prop.cultg / (prop.crop + prop.cultg))^2) * (ocg.sd^2)
) * g_scalar_oc


# calculate initial bd value
# note that g_scalar_bd is applied within the calc.bd function
bd <- calc.bd(exrisk.t0.slope)
#bd <- calc.bd(exrisk.t0)   # <------------ testing purposes only
summary(bd)



# ### SET RESTORATION AREA TARGET ###
# 
# # we should set the total area target once only and use this for all scenarios
# 
# # ******* note there are multiple, alternative targets here - just set one of them:  ********
# 
# # 350 million hectares (NY Declaration target)  "t350" #########################
# restoration.area <- 3.5E6
# # do not change this line:
# restoration.rhs <- (restoration.area / A) * g_scalar_area
# 
# # set the output directory for to which all derived data will be written or overwritten
# # must end in /
# outdir <- "opt_results_t350_v5/"
# if (!dir.exists(outdir)) dir.create(outdir)
# 
# 
# # 150 million hectares (NY Declaration target)  "t150" ######################### 
# #restoration.area <- 1.5E6
# # do not change this line:
# restoration.rhs <- (restoration.area / A) * g_scalar_area
# 
# # set the output directory for to which all derived data will be written or overwritten
# # must end in /
# #outdir <- "opt_results_t150_v5/"
# if (!dir.exists(outdir)) dir.create(outdir)

# need to load this for the plotting functions:
load(paste0(dir, "terrestrial_lands.RData"))

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg) * ub.perc.constraint



### PARAMETER DEFINITIONS
# number of cores to use in parallel processing
ncores <- 4


### CAUTION: important to reset oc after running the postprocessing code, so we just reset it
# here to be safe:
oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg 
oc <- oc * g_scalar_oc


### run all the scenarios:
#source("run_scenarios.r")
source("wrap_scenarios.r")

#if (g_fatal_error) message("FATAL ERROR DETECTED")




# ### SCENARIO: WHOLE WORLD (MAX RESTORATION) ####################################################
# 
# # this is not a target-specific scenario, so it only needs to be run a single time among all targets
# 
# scen <- "scen_world"
# w <- "NA"
# nsteps <- 1
# 
# # calculate initial bd value - we reset this every time as a precautionary measure
# bd <- calc.bd(exrisk.t0.slope)
# 
# result <- list()
# result$x <- ub
# 
# # create the same data structures as are output from the optimisations
# res.prop.restored.pu <- matrix(0, nrow=np, ncol=nsteps)
# res.total.restored.pu <- rep(0, np)
# res.area.restored.spp <- matrix(0, nrow=ns, ncol=nsteps)
# res.total.restored.spp <- rep(0, ns)
# res.exrisk <- matrix(0, nrow=ns, ncol=nsteps)
# res.objval <- rep(0, nsteps)
# delta.hab.pu <- rep(0, np) 
# 
# delta.hab.pu <- result$x
# # allocate that restoration to habitat types
# delta.vegtype.pu <- prop.restore * delta.hab.pu
# 
# # convert the change in habitat per pu to change in habitat per species
# delta.hab.spp <- rep(0, ns)
# for (i in 1:ns){
# 	for(j in 1:5){
# 		if (sphabm_proc[i,j] == 1){
# 			delta.hab.spp[i] <- delta.hab.spp[i] + sum(delta.vegtype.pu[species_index_list_proc[[i]], j], na.rm=T)
# 		}
# 	}
# }
# delta.hab.spp <- delta.hab.spp * A
# res.area.restored.spp[,1] <- delta.hab.spp
# res.total.restored.spp <- res.total.restored.spp + delta.hab.spp
# 
# res.exrisk[,1] <- extinction.risk(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
# res.prop.restored.pu[,1] <- result$x
# res.total.restored.pu <- result$x
# 
# save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
# save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
# save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
# save(res.exrisk, file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))




















# 
# ### NOTES AND TEMPORARY CODE: IGNORE EVERYTHING BELOW THIS POINT! #########################
# 
# 
# ### DIAGNOSTICS TO LOOK AT BD PLATEAU PROBLEM ###################
# 
# # these BD indices have the same max BD value at t0:
# q <- c(214524, 3159507)
# # corresponding to these 
# qq <- master_index[q]
# # cell ids: 2918594 20964361
# 
# q.sp1 <- rep(0, ns)
# q.sp2 <- rep(0, ns)
# 
# for (i in 1:ns){
# 	if (!is.null(species_index_list[[i]])){
# 		if (length(which(species_index_list[[i]] == q[1])) == 1) q.sp1[i] <- 1
# 		if (length(which(species_index_list[[i]] == q[2])) == 1) q.sp2[i] <- 1
# 	}
# }
# sum(q.sp1)
# sum(q.sp2)
# 
# 
# 
# # bd value (starting)
# bda <- rep(0, np)
# for (i in 1:ns){
# 	bda <- bda + exrisk.t0.slope[i] * (prop.restore %*% sphabm_proc[i,])
# }
# bda <- bda * g_scalar_bd
# summary(bda)
# sum(bda - bd.t0)
# cor(bda - bd.t0)
# 
# 
# # bd value (starting)
# Sys.time()
# bda2 <- rep(0, np)
# for (i in 1:ns){
# 	if (!is.null(species_index_list[[i]])){
# 		bda2[species_index_list[[i]]] <- bda2[species_index_list[[i]]] + exrisk.t0.slope[i] * (prop.restore[species_index_list[[i]],] %*% sphabm_proc[i,])
# 	}
# }
# bda2 <- bda2 * g_scalar_bd
# Sys.time()
# 
# 
# summary(bda2)
# sum(bda2 - bd.t0)
# cor(bda2 - bd.t0)
# 
# 
# # plot map of objective function
# plot.pu.map(bda2, fname=paste0(outdir, scen, "_map.fixtest_1.png"))
# 
# 
# 
# # calculate initial bd value
# Sys.time()
# bd <- rep(0, np)
# for (i in 1:dim(usphab_proc)[1]){
# 	hab.values <- prop.restore %*% usphab_proc[i,]
# 	for (j in 1:length(usphab_index[[i]])){
# 		if (!is.null(species_index_list[[usphab_index[[i]][j]]])){
# 			recs <- species_index_list[[usphab_index[[i]][j]]]
# 			bd[recs] <- bd[recs] + exrisk.t0.slope[usphab_index[[i]][j]] * hab.values[recs]
# 		}
# 	}
# 	
# }
# bd <- bd * g_scalar_bd
# Sys.time()
# summary(bd)
# 
# cnt <- 0
# for (j in 1:length(species_index_list)){
# 	if (is.null(species_index_list[[j]])) cnt <- cnt + 1
# }
# cnt
# 
# # make OC a raster for Bernardo:
# 
# r <- r.er
# values(r) <- rep(NA, ncell(r))
# values(r)[master_index] <- oc
# writeRaster(r, file="oc.tif")
