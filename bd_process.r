# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

### RUN AUXILIARY ANALYSES (ORIGINAL AREA AND OPPORTUNITY COSTS)
aux.run = F

if (aux.run){
  source('~/Documents/IIS_PROJECTS/global_rest_prior/OA/OA_preprocessing.R')
  source('~/Documents/IIS_PROJECTS/global_rest_prior/OA/OA_calc.R')
  rm(list=ls())
  
  source('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_econ/econ-definitions.R')
  source('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_econ/econ-oc.R')
  rm(list=ls())  
}


# load libraries
library(sp)
library(raster)

# set folder into which all pre-processed datasets will be written or
# overwritten if they already exist, must end with a /, and must already exist
# on disk
dir = 'inputdata_v5/'
setwd('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_priorization/rawdata/')
source("../functions.r")


### LOAD BASE rcomposeS

# (ZERO) terrestrial lands
r.wat = raster("./900m/current_LU/ESA_landuse_300m_2015_water_media_900m_Molweide.tif")
r.terr = abs(r.wat - 1) > 1.e-7

A <- ((res(r.terr)[1])^2) / 1E6
save(A, file=paste0(dir, "A_900.RData"))


# Land-use-based rasters are not aligned with species rasters in 900m resolution
# we need to correct for that. We use the function rcompose to create rasters
# that are already intersected to the correct extent, r.ext, as computed below
r.ext = extent(intersect(r.terr, rcompose('./900m/BD_ranges/amphibians/100099193.tif')))


# (A) proportions of landcovers (current)
r.plc1 <- rcompose("./900m/current_LU/crop_class11_2015_900m_Moll.tif")
r.plc2 <- rcompose("./900m/current_LU/CultGrass_class11_2015_900m_Moll.tif")
r.plc3 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_forest_media_900m_Molweide.tif")
r.plc4 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_wetlands_media_900m_Molweide.tif")
r.plc5 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_desert_media_900m_Molweide.tif")
r.plc6 <- rcompose("./900m/current_LU/NatGrass_2015_900m_Moll.tif")
r.plc7 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_shrubland_media_900m_Molweide.tif")
r.plc8 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_ice_media_900m_Molweide.tif") + rcompose(r.wat)
#r.plc11 <- rcompose("./900m/current_LU/ESA_landuse_300m_2015_Moll_class11_c_byte.tif")

# correcting ESA's class-11 reclassification to half-pasture / half-crop mixture
#r.plc1 = r.plc1 - (0.5*c11)
#r.plc2 = r.plc2 + (0.5*c11)

# (B) proportions of landcovers (original)
# source('./900m/original_LC/renormOA.R')
r.oa1 <- rcompose("./900m/original_LC/OAforest900m_renorm.tif")
r.oa2 <- rcompose("./900m/original_LC/OAwetland900m_renorm.tif")
r.oa3 <- rcompose("./900m/original_LC/OAdesert900m_renorm.tif")
r.oa4 <- rcompose("./900m/original_LC/OAgrassland900m_renorm.tif")
r.oa5 <- rcompose("./900m/original_LC/OAshrubland900m_renorm.tif")
r.oa8 <- rcompose("./900m/original_LC/OAicewaterrock900m_renorm.tif")


# (C) elevation
r.elev <- rcompose("./900m/DEM/DEM_World_900m_Molweide3.tif")


# (D) species range rcomposes
spp_dir <- "./900m/BD_ranges/mammals/"
rnames.m <- list.files(spp_dir, full.names=TRUE, pattern = "\\.tif$")
length(rnames.m)

spp_dir <- "./900m/BD_ranges/amphibians"
rnames.a <- list.files(spp_dir, full.names=TRUE, pattern = "\\.tif$")
length(rnames.a)

spp_dir <- "./900m/BD_ranges/birds"
rnames.b <- list.files(spp_dir, full.names=TRUE, pattern = "\\.tif$")
length(rnames.b)

spp_raster_names <- c(rnames.m, rnames.b, rnames.a)
spp_taxon <- c(rep("M", length(rnames.m)), rep("A", length(rnames.a)), rep("B", length(rnames.b)))
length(spp_raster_names)

# extract species numbers
spid <- rep(0, length(spp_raster_names))
for (i in 1:length(spp_raster_names)){
  s <- strsplit(spp_raster_names[i], "/")[[1]]
  s2 <- s[length(s)]
  spid[i] <- as.numeric(substr(s2, 1, nchar(s2) - 4))
}

# check for NAs:
summary(spid)

# error check for duplicate species IDs
if(!(length(spid) == length(unique(spid)))) message("ERROR: duplicate IDs")





# create maps of the summed species ranges for each taxon and the total - useful for error checking
# v.m <- rep(0, ncell(r.er))
# v.a <- rep(0, ncell(r.er))
# v.b <- rep(0, ncell(r.er))
# 
# for (i in 1:length(spp_raster_names)){
# 	if (spp_taxon[i] == "M"){
# 		v.m <- v.m + values(rcompose(spp_raster_names[i]))
# 	} else {
# 		if (spp_taxon[i] == "A"){
# 			v.a <- v.a + values(rcompose(spp_raster_names[i]))
# 		} else {
# 			v.b <- v.b + values(rcompose(spp_raster_names[i]))
# 		}
# 	}
# }
# 
# v.terr <- values(r.terr)
# r <- rcompose(spp_raster_names[1])
# values(r) <- replace(v.m, which(v.m == 0 | v.terr == 0 | is.na(v.terr)), NA)
# writercompose(r, file=paste0(dir, "spp_sum_mammals.tif"), overwrite=T)
# colfunc <- colorRampPalette(c("green4", "green", "yellow", "orange", "red", "red4"), bias=1)
# cols <- colfunc(32)
# png(file="map_spp_sum_mammals.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, col=c("white", "grey90"), axes=F, box=F, legend=FALSE)
# plot(r, col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# values(r) <- replace(v.a, which(v.a == 0 | v.terr == 0 | is.na(v.terr)), NA)
# writercompose(r, file=paste0(dir, "spp_sum_amphibians.tif"), overwrite=T)
# png(file="map_spp_sum_amphibians.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, col=c("white", "grey90"), axes=F, box=F, legend=FALSE)
# plot(r, col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# values(r) <- replace(v.b, which(v.b == 0 | v.terr == 0 | is.na(v.terr)), NA)
# writercompose(r, file=paste0(dir, "spp_sum_birds.tif"), overwrite=T)
# png(file="map_spp_sum_birds.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, col=c("white", "grey90"), axes=F, box=F, legend=FALSE)
# plot(r, col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# v.all <- v.m + v.a + v.b
# values(r) <- replace(v.all, which(v.all == 0 | v.terr == 0 | is.na(v.terr)), NA)
# writercompose(r, file=paste0(dir, "spp_sum_all.tif"), overwrite=T)
# png(file="map_spp_sum_all.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.terr, col=c("white", "grey90"), axes=F, box=F, legend=FALSE)
# plot(r, col=cols, axes=F, box=F, add=T)
# dev.off()
# 
# rm(v.m)
# rm(v.a)
# rm(v.b)
# rm(v.all)



### MASTER LIST OF CELL INDICES
# create master list of cell indices from ecoregion dataset (or other suitable dataset)

# master index for all terrestrial cells (except rock) used for the purpose of calculating total, global species ranges
#terrestrial_index <- which(!(is.na(values(r.er)) | values(r.er) == 0 | values(r.er) == 11))
terrestrial_index = which(values(r.terr) > 1.e-6)
length(terrestrial_index)
save(terrestrial_index, file=paste0(dir, "terrestrial_index_900.RData"))

# master index for cells that are part of the optimisation, after removing cells with no restorable area, and any that have NoData values in the reference datasets
#master_index <- which(!(is.na(values(r.er)) | values(r.er) == 0 | values(r.er) == 11))
master_index = terrestrial_index
save(master_index, file=paste0(dir, "master_index_900.RData"))
length(master_index)
# 5037306

load(file=paste0(dir, "terrestrial_index_900.RData"))
load(file=paste0(dir, "master_index_900.RData"))



### VECTOR OF UPPER BOUND OF DECISION VARIABLES (MAX OF EACH CELL THAT CAN BE RESTORED)

# there were formerly problems with NoData and negative values in these rcomposes, so we need to trap for that
# (this code does nothing if there are no NA values - so it is safe to run)
prop.crop <- values(r.plc1)[master_index]
recs <- which(is.na(prop.crop))
length(recs)
if (length(recs) > 0) prop.crop[recs] <- 0

recs <- which(prop.crop < 0)
length(recs)
if (length(recs) > 0) prop.crop[recs] <- 0

prop.cultg <- values(r.plc2)[master_index]
recs <- which(is.na(prop.cultg))
length(recs)

recs <- which(prop.cultg < 0)
length(recs)
if (length(recs) > 0) prop.cultg <- replace(prop.cultg, which(prop.cultg < 0), 0)


# calculate upper bound of decision variables
upperbound <- prop.crop + prop.cultg
summary(upperbound)

# we need to remove cells that cannot be restored to reduce the size of the problem, 
# so we remove those cells from master_index
#recs <- which(upperbound < 0.01)
recs <- which(upperbound == 0)
length(recs)
if (length(recs) > 0){
  upperbound <- upperbound[-recs]
  master_index <- master_index[-recs]
}

min(upperbound)
length(master_index)



# remove any cells that have NAs for any other dataset
rid <- c()
vals <- values(r.plc1)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc2)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc3)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc4)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc5)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc6)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc7)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.plc8)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
# vals <- values(r.plc9)[master_index]
# rid <- c(rid, which(is.na(vals)))
# length(rid)

vals <- values(r.oa1)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.oa2)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.oa3)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.oa4)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)
vals <- values(r.oa5)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)

vals <- values(r.elev)[master_index]
rid <- c(rid, which(is.na(vals)))
length(rid)

rid <- unique(rid)
master_index <- master_index[-rid]
length(master_index)

save(master_index, file=paste0(dir, "master_index_900.RData"))


load(file=paste0(dir, "master_index_900.RData"))


# now recalculate upperbound with the new master index
prop.crop <- values(r.plc1)[master_index]
recs <- which(is.na(prop.crop))
length(recs)
if (length(recs) > 0) prop.crop[recs] <- 0

recs <- which(prop.crop < 0)
length(recs)
if (length(recs) > 0) prop.crop[recs] <- 0

prop.cultg <- values(r.plc2)[master_index]
recs <- which(is.na(prop.cultg))
length(recs)

recs <- which(prop.cultg < 0)
length(recs)
if (length(recs) > 0) prop.cultg <- replace(prop.cultg, which(prop.cultg < 0), 0)


# calculate upper bound of decision variables
upperbound <- prop.crop + prop.cultg
summary(upperbound)
min(upperbound)

save(upperbound, file=paste0(dir, "upperbound.RData"))



### ERROR CHECKING
# confirm all input rcomposes have the same dimensions and projection, and
# projection is defined this takes a few minutes to run for the species data, so
# no need to run this if you have already run this before and the data has not
# changed

# check.params <- function(r, r.master){
#   if(!(nrow(r) == nrow(r.master) & ncol(r) == ncol(r.master) & projection(r) == projection(r.master))){
#     return(1)
#   } else {
#     return(0)
#   }
# }
# 
# if (check.params(r.elev, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc1, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc2, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc3, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc4, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc5, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc6, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc7, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.plc8, r.terr) == 1) message("ERROR: mismatch")
# #if (check.params(r.plc9, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.oa1, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.oa2, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.oa3, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.oa4, r.terr) == 1) message("ERROR: mismatch")
# if (check.params(r.oa5, r.terr) == 1) message("ERROR: mismatch")
# 
# for (i in 1:length(spp_raster_names)){
#   if (i %% 500 == 0) message(i)
#   r.sp <- rcompose(spp_raster_names[i])
#   if (check.params(r.sp, r.terr) == 1){
#     message(paste0("ERROR: ", spp_raster_names[i], " ", nrow(r.sp), " ", ncol(r.sp)))
#   }
# }

# check proportions sum to 1: current land cover
vals <- values(r.plc1)[master_index]
vals <- vals + values(r.plc2)[master_index]
vals <- vals + values(r.plc3)[master_index]
vals <- vals + values(r.plc4)[master_index]
vals <- vals + values(r.plc5)[master_index]
vals <- vals + values(r.plc6)[master_index]
vals <- vals + values(r.plc7)[master_index]
vals <- vals + values(r.plc8)[master_index]
#vals <- vals + values(r.plc9)[master_index]
summary(vals)
length(which(!vals == 1))
length(which(!vals == 1))/length(master_index)

# check proportions sum to 1: current land cover
vals <- values(r.oa1)[master_index]
vals <- vals + values(r.oa2)[master_index]
vals <- vals + values(r.oa3)[master_index]
vals <- vals + values(r.oa4)[master_index]
vals <- vals + values(r.oa5)[master_index]
summary(vals)
length(which(!vals == 1))
length(which(!vals == 1))/length(master_index)



### ELEVATION VALUE

elev <- values(r.elev)[master_index]
summary(elev)
save(elev, file=paste0(dir, "elev.RData"))
summary(elev)


# note that because elevation is used to calculate the current and previous
# potential range extents, we need one elevation vector associated with the
# master_index (elev) and one associated with the terrestrial_index:
elev.terr <- values(r.elev)[terrestrial_index]
summary(elev.terr)
# note that we cannot replace the NA's here, so we just have to ignore them in the potential species range calculations
save(elev.terr, file=paste0(dir, "elev.terr.RData"))
summary(elev.terr)


# PROP OF CROP AND CULT GRASS

prop.crop <- values(r.plc1)[master_index]
prop.cultg <- values(r.plc2)[master_index]
save(prop.crop, file=paste0(dir, "prop.crop.RData"))
save(prop.cultg, file=paste0(dir, "prop.cultg.RData"))


###	LIST OF CELL INDICES FOR EACH SPECIES - WARNING - SLOW AND LARGE
# for each species, record the cell indices that are non-zero and that are within the set of master indices (eliminates all ocean cells)
# This MUST refer to the index of the master_index records, not the cell indices of the range rcomposes directly!

# This is a memory efficient but complex data structure and it is easy to break
# this. The master_index is an index of cell positions in the rcomposes, and the
# species_index_list is a list of indices within the master_index (hence an
# index of an index). It is essential to take care to ensure that the right
# index is associated with the right data structures.

# load elevation limits data - but note that they are semicolon delimited
elim.lower <- read.csv("./5km/others/TB_LowerElevation_MAB_May2018.csv", sep=";")
elim.upper <- read.csv("./5km/others/TB_UpperElevation_MAB_May2018.csv", sep=";")
head(elim.lower)
head(elim.upper)

# optional:
# do some quick checks to determine what species have elevations
has.elev.lower <- rep(0, length(spid))
has.elev.upper <- rep(0, length(spid))
for (i in 1:length(spid)){
  rec <- which(elim.lower$taxonid == spid[i])
  if (length(rec) == 1){
    if (!is.na(elim.lower$lowerelevation[rec])) has.elev.lower[i] <- 1
  } 
  rec <- which(elim.upper$taxonid == spid[i])
  if (length(rec) == 1){
    if (!is.na(elim.upper$upperelevation[rec])) has.elev.upper[i] <- 1
  } 
}
sum(has.elev.lower)
sum(has.elev.upper)
table(spp_taxon, has.elev.lower)
table(spp_taxon, has.elev.upper)


load(file=paste0(dir, "master_index.RData"))

species_index_list <- list()
errors <- rep(0, length(spp_raster_names))
species_cell_cnt <- rep(0, length(spp_raster_names))

for (i in 1:length(spp_raster_names)){
  if (i %% 1000 == 0){
    message(i)
    save(species_index_list, file=paste0(dir, "species_index_list.RData"))
    save(species_cell_cnt, file=paste0(dir, "species_cell_cnt.RData"))
  }
  #if (class(try(values(rcompose(spp_raster_names[i])), silent=T)) == "try-error"){
  if (class(try(rcompose(spp_raster_names[i]), silent=T)) == "try-error"){
    spp_raster_names = spp_raster_names[-i] # removes species from list of rcompose names
    spid = spid[-i]
    print(paste(spid[i], 'removed'))
    next                                    # jumps to the next i-th iteration
  }
  r.sp <- rcompose(spp_raster_names[i])
  v <- values(r.sp)[master_index]
  rec.l <- which(elim.lower$taxonid == spid[i])
  rec.u <- which(elim.upper$taxonid == spid[i])
  elev.upperlim <- NA
  elev.lowerlim <- NA
  if (length(rec.u) == 1) elev.upperlim <- elim.upper$upperelevation[rec.u]
  if (length(rec.l) == 1) elev.lowerlim <- elim.lower$lowerelevation[rec.l]
  
  if (is.na(elev.upperlim) & is.na(elev.lowerlim)){
    recs <- which(v == 1)
  } else {
    if (!is.na(elev.upperlim) & !is.na(elev.lowerlim)){
      recs <- which(v == 1 & elev <= elev.upperlim & elev >= elev.lowerlim)
    } else {
      if (!is.na(elev.upperlim)){
        recs <- which(v == 1 & elev <= elev.upperlim)
      } else {
        recs <- which(v == 1 & elev >= elev.lowerlim)
      }
    }
  }
  
  species_cell_cnt[i] <- length(recs)
  if (length(recs) == 0){
    species_index_list[[i]] <- NULL
    errors[i] <- i
  } else {
    species_index_list[[i]] <- recs
  }
} # i

save(species_index_list, file=paste0(dir, "species_index_list.RData"))
save(species_cell_cnt, file=paste0(dir, "species_cell_cnt.RData"))




# repeat this for all terrestrial cells for purpose of calculating t0 baseline

load(file=paste0(dir, "terrestrial_index.RData"))
load(file=paste0(dir, "elev.terr.RData"))

species_index_list_terr <- list()
errors <- rep(0, length(spp_raster_names))
species_cell_cnt_terr <- rep(0, length(spp_raster_names))

for (i in 1:length(spp_raster_names)){
  if (i %% 1000 == 0){
    message(i)
    save(species_index_list_terr, file=paste0(dir, "species_index_list_terr.RData"))
    save(species_cell_cnt_terr, file=paste0(dir, "species_cell_cnt_terr.RData"))
  }
  r.sp <- rcompose(spp_raster_names[i])
  v <- values(r.sp)[terrestrial_index]
  rec.l <- which(elim.lower$taxonid == spid[i])
  rec.u <- which(elim.upper$taxonid == spid[i])
  elev.upperlim <- NA
  elev.lowerlim <- NA
  if (length(rec.u) == 1) elev.upperlim <- elim.upper$upperelevation[rec.u]
  if (length(rec.l) == 1) elev.lowerlim <- elim.lower$lowerelevation[rec.l]
  
  if (is.na(elev.upperlim) & is.na(elev.lowerlim)){
    recs <- which(v == 1)
  } else {
    if (!is.na(elev.upperlim) & !is.na(elev.lowerlim)){
      recs <- which(v == 1 & elev.terr <= elev.upperlim & elev.terr >= elev.lowerlim)
    } else {
      if (!is.na(elev.upperlim)){
        recs <- which(v == 1 & elev.terr <= elev.upperlim)
      } else {
        recs <- which(v == 1 & elev.terr >= elev.lowerlim)
      }
    }
  }
  
  species_cell_cnt_terr[i] <- length(recs)
  if (length(recs) == 0){
    species_index_list_terr[[i]] <- NULL
    errors[i] <- i
  } else {
    species_index_list_terr[[i]] <- recs
  }
}

save(species_index_list_terr, file=paste0(dir, "species_index_list_terr.RData"))
save(species_cell_cnt_terr, file=paste0(dir, "species_cell_cnt_terr.RData"))



load(file=paste0(dir, "species_index_list_terr.RData"))
load(file=paste0(dir, "species_cell_cnt_terr.RData"))



# check the counts


length(which(species_cell_cnt == 0))
length(which(species_cell_cnt_terr == 0))
length(species_cell_cnt_terr)



### ADJUST THE SPECIES INDEX LIST TO REMOVE ADDITIONAL CELLS

# only run this if the master_index has changed prior to initially creating the species_index_list

# species_index_list_rev <- list()

# for (i in 1:length(species_index_list)){
# 	if (is.null(species_index_list[[i]])){
# 		species_index_list_rev[[i]] <- NULL
# 	} else {
# 		recs <- which(species_index_list[[i]] %in% master_index)
# 		if (length(recs) > 0){
# 			species_index_list_rev[[i]] <- species_index_list[[i]][recs]
# 		} else {
# 			species_index_list_rev[[i]] <- NULL
# 		}
# 	}
# }

# save(species_index_list_rev, file=paste0(dir, "species_index_list_rev.RData"))

# object.size(species_index_list_rev)


load(file=paste0(dir, "species_index_list.RData"))






### SPECIES X HABITAT MATRIX AND CONSERVATION STATUS

sphabdf <- read.csv("./rawdata/5km/BD_ranges/Habitats_marginal_excluded_habitats_selected.csv")
head(sphabdf)

unique(sphabdf$habitat)
sphabdf[which(sphabdf$habitat == 3)[1],]


# codes:
# 1: forest
# 2: savannah
# 3: shruband
# 4: grassland
# 5: wetlands
# 6: rocky areas
# 8: desert
# 14: Artificial/Terrestrial - Subtropical/Tropical Heavily Degraded Former Forest suitability season majorimportance

# codes for species x habitat matrix (numbers refer to column indices:
# 1: forest
# 2: wetland
# 3: desert
# 4: natural grassland
# 5: shrubland

# conservation status:
# > unique(sphabdf$category)
# [1] "LC" "NT" "DD" "VU" "EN" "CR" "EW" "EX"

# 1: CR
# 2: EN
# 3: NT
# 4: VU
# 5: LC
# 10: EX
# 11: EW
# 12: DD
# 13: NE 

# Extinct (EX) – No known individuals remaining
# Extinct in the wild (EW) – Known only to survive in captivity, or as a naturalized population outside its historic range
# Critically endangered (CR) – Extremely high risk of extinction in the wild
# Endangered (EN) – High risk of extinction in the wild
# Vulnerable (VU) – High risk of endangerment in the wild
# Near threatened (NT) – Likely to become endangered in the near future
# Least concern (LC) – Lowest risk (Does not qualify for a more at-risk category; widespread and abundant taxa are included in this category.)
# Data deficient (DD) – Not enough data to make an assessment of its risk of extinction
# Not evaluated (NE) – Has not yet been evaluated against the criteria






# populate species x habitat matrix

sphabm <- matrix(0, nrow=length(spid), ncol=5)
spcat <- rep(0, length(spid))
errors <- c()
sphabclass <- rep("M", dim(sphabdf)[1])
sphabclass[which(sphabdf$class_name == "AVES")] <- "B"
sphabclass[which(sphabdf$class_name == "AMPHIBIA")] <- "A"

# just test all spid are unique
(length(spid) == length(unique(spid)))

for (i in 1:length(spid)){
  #recs <- which(sphabdf$taxonid == spid[i] & sphabclass == spp_taxon[i])
  recs <- which(sphabdf$taxonid == spid[i])
  
  if (length(recs) == 0){
    #message(paste0("No records found for species ID: ", spid[i]))
    if (!spp_taxon[i] == "R") errors <- c(errors, spid[i])
  } else {
    # encode species category
    if (sphabdf$category[recs[1]] == "CR") spcat[i] <- 1
    if (sphabdf$category[recs[1]] == "EN") spcat[i] <- 2
    if (sphabdf$category[recs[1]] == "NT") spcat[i] <- 3
    if (sphabdf$category[recs[1]] == "VU") spcat[i] <- 4
    if (sphabdf$category[recs[1]] == "LC") spcat[i] <- 5
    if (sphabdf$category[recs[1]] == "EX") spcat[i] <- 10
    if (sphabdf$category[recs[1]] == "EW") spcat[i] <- 11
    if (sphabdf$category[recs[1]] == "DD") spcat[i] <- 12
    if (sphabdf$category[recs[1]] == "NE") spcat[i] <- 13
    
    # check for forest suitability
    rec <- which(sphabdf$habitat[recs] == 1)
    if (length(rec) > 0) sphabm[i,1] <- 1
    # check for wetland suitability
    rec <- which(sphabdf$habitat[recs] == 5)
    if (length(rec) > 0) sphabm[i,2] <- 1
    # check for desert suitability
    rec <- which(sphabdf$habitat[recs] == 8)
    if (length(rec) > 0) sphabm[i,3] <- 1
    # check for natural grassland suitability
    rec <- which(sphabdf$habitat[recs] == 4)
    if (length(rec) > 0) sphabm[i,4] <- 1
    # check for shrubland suitability
    rec <- which(sphabdf$habitat[recs] == 2 | sphabdf$habitat[recs] == 3)
    if (length(rec) > 0) sphabm[i,5] <- 1
  }
}

dim(sphabm)
length(errors)
save(errors, file=paste0(dir, "spid_no_habitat_lookup_entires.RData"))


# check that every species has at least 1 suitable habitat
errors2 <- c()
for (i in 1:dim(sphabm)[1]){
  if (sum(sphabm[i,]) == 0) errors2 <- c(errors2, spid[i])
}
length(errors2)
save(errors2, file=paste0(dir, "spid_no_habitat_lookup_entires2.RData"))

# check for impermissable categories: extinct, extinct in wild, etc
length(which(spcat > 5))

recs <- which(spid %in% errors2)
table(spcat[recs])

recs <- which(species_cell_cnt == 0)
table(spcat[recs])





### CREATE LIST OF SPECIES INDICES FOR PROCESSING 

# rather than subset all the data to remove problematic species, which is hard because there are so many interrelated data structures, we just maintain a list of indices representing species that will be used in the optimisation
# this is also beneficial because different scenarios use different sets of species


valid_species <- rep(TRUE, dim(sphabm)[1])

# remove species with no cells
for (i in 1:length(species_index_list)){
  if (is.null(species_index_list[[i]])) valid_species[i] <- FALSE
}

# remove species with sphabm errors
for (i in 1:dim(sphabm)[1]){
  if (sum(sphabm[i,]) == 0) valid_species[i] <- FALSE
}


# we need to identify the subset of species that have no presettlement habitat
# what about the 600+ species that have zero pre-settlement habitat?

prop.hab.oa <- matrix(0, nrow=length(terrestrial_index), ncol=5)
v <- values(r.oa1)[terrestrial_index]
prop.hab.oa[,1] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa2)[terrestrial_index]
prop.hab.oa[,2] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa3)[terrestrial_index]
prop.hab.oa[,3] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa4)[terrestrial_index]
prop.hab.oa[,4] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa5)[terrestrial_index]
prop.hab.oa[,5] <- replace(v, which(is.na(v)), 0)

habarea.max <- rep(0, dim(sphabm)[1])
for (i in 1:dim(sphabm)[1]){
  for(j in 1:5){
    if (sphabm[i,j] == 1){
      habarea.max[i] <- habarea.max[i] + sum(prop.hab.oa[species_index_list_terr[[i]], j], na.rm=T)
    }
  }
}
summary(habarea.max)

recs <- which(habarea.max == 0)
length(recs)
valid_species[recs] <- FALSE

# destroy these data structures as they are not correct and are recreated later:
rm(habarea.max)
rm(prop.hab.oa)

# remove extinct and extinct in wild
recs <- which(spcat == 10 | spcat == 11)
length(recs)
valid_species[recs] <- FALSE

table(valid_species)

table(valid_species, spcat)
#table(valid_species, spp_taxon)





# reduce species data structures to smallest set required for processing
recs <- which(valid_species)
length(recs)
sphabm_proc <- sphabm[recs,]
spid_proc <- spid[recs]
spcat_proc <- spcat[recs]
spp_taxon_proc <- spp_taxon[recs]

save(sphabm_proc, file=paste0(dir, "sphabm_proc.RData"))
save(spid_proc, file=paste0(dir, "spid_proc.RData"))
save(spcat_proc, file=paste0(dir, "spcat_proc.RData"))
save(spp_taxon_proc, file=paste0(dir, "spp_taxon_proc.RData"))


# subset the species index data structures to remove species we cannot process
load(file=paste0(dir, "species_index_list_terr.RData"))

species_index_list_terr_proc <- list()
for (i in 1:length(recs)){
  species_index_list_terr_proc[[i]] <- species_index_list_terr[[recs[[i]]]]
}
save(species_index_list_terr_proc, file=paste0(dir, "species_index_list_terr_proc.RData"))


load(file=paste0(dir, "species_index_list.RData"))

species_index_list_proc <- list()
for (i in 1:length(recs)){
  species_index_list_proc[[i]] <- species_index_list[[recs[[i]]]]
}
save(species_index_list_proc, file=paste0(dir, "species_index_list_proc.RData"))





### CALCULATE Amax, THE MAXIMUM POTENTIAL HABITAT (used in the extinction risk calculation)

load(file=paste0(dir, "terrestrial_index.RData"))
load(file=paste0(dir, "species_index_list_proc.RData"))
load(file=paste0(dir, "species_index_list_terr_proc.RData"))
load(file=paste0(dir, "sphabm_proc.RData"))
load(file=paste0(dir, "spid_proc.RData"))
load(file=paste0(dir, "spcat_proc.RData"))

prop.hab.oa <- matrix(0, nrow=length(terrestrial_index), ncol=5)
v <- values(r.oa1)[terrestrial_index]
prop.hab.oa[,1] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa2)[terrestrial_index]
prop.hab.oa[,2] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa3)[terrestrial_index]
prop.hab.oa[,3] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa4)[terrestrial_index]
prop.hab.oa[,4] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa5)[terrestrial_index]
prop.hab.oa[,5] <- replace(v, which(is.na(v)), 0)

habarea.max <- rep(0, dim(sphabm_proc)[1])
for (i in 1:dim(sphabm_proc)[1]){
  for(j in 1:5){
    if (sphabm_proc[i,j] == 1){
      habarea.max[i] <- habarea.max[i] + sum(prop.hab.oa[species_index_list_terr_proc[[i]], j], na.rm=T)
    }
  }
}
habarea.max <- habarea.max * A
summary(habarea.max)
min(habarea.max)


### CALCULATE AT0 (CURRENT HABITAT AREA)

# At0
prop.hab.t0 <- matrix(0, nrow=length(terrestrial_index), ncol=5)
v <- values(r.plc3)[terrestrial_index]
prop.hab.t0[,1] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc4)[terrestrial_index]
prop.hab.t0[,2] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc5)[terrestrial_index]
prop.hab.t0[,3] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc6)[terrestrial_index]
prop.hab.t0[,4] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc7)[terrestrial_index]
prop.hab.t0[,5] <- replace(v, which(is.na(v)), 0)

# area of habitat for each species = habitat suitability * habitat availability * range
habarea.t0 <- rep(0, dim(sphabm_proc)[1])
for (i in 1:dim(sphabm_proc)[1]){
  for(j in 1:5){
    if (sphabm_proc[i,j] == 1){
      habarea.t0[i] <- habarea.t0[i] + sum(prop.hab.t0[species_index_list_terr_proc[[i]], j], na.rm=T)
    }
  }
}
habarea.t0 <- habarea.t0 * A
summary(habarea.t0)
min(habarea.t0)


summary(habarea.t0)
summary(habarea.max)
summary(prop.hab.t0)
#summary(prop.hab.max)

# extinction risk
exrisk.t0 <- extinction.risk(habarea.t0, habarea.max, z=0.25)
summary(exrisk.t0)

# extinction risk slopes
exrisk.t0.slope <- calc.extinction.slope(habarea.t0, habarea.max, z=0.25)
summary(exrisk.t0.slope)


# why large slopes?
recs <- which(exrisk.t0.slope > 1)
habarea.max[recs]
habarea.t0[recs]


# proportion hab within master index cells only
prop.hab.pu <- matrix(0, nrow=length(master_index), ncol=5)
v <- values(r.plc3)[master_index]
prop.hab.pu[,1] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc4)[master_index]
prop.hab.pu[,2] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc5)[master_index]
prop.hab.pu[,3] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc6)[master_index]
prop.hab.pu[,4] <- replace(v, which(is.na(v)), 0)
v <- values(r.plc7)[master_index]
prop.hab.pu[,5] <- replace(v, which(is.na(v)), 0)


save(prop.hab.pu, file=paste0(dir, "prop.hab.pu.RData"))
save(prop.hab.t0, file=paste0(dir, "prop.hab.t0.RData"))
save(prop.hab.oa, file=paste0(dir, "prop.hab.oa.RData"))
save(habarea.t0, file=paste0(dir, "habarea.t0.RData"))
save(habarea.max, file=paste0(dir, "habarea.max.RData"))
save(exrisk.t0, file=paste0(dir, "exrisk.t0.RData"))
save(exrisk.t0.slope, file=paste0(dir, "exrisk.t0.slope.RData"))

# when restoration occurs it happens in proportion to the frequencies in OA

prop.restore <- matrix(0, nrow=length(master_index), ncol=5)
v <- values(r.oa1)[master_index]
prop.restore[,1] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa2)[master_index]
prop.restore[,2] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa3)[master_index]
prop.restore[,3] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa4)[master_index]
prop.restore[,4] <- replace(v, which(is.na(v)), 0)
v <- values(r.oa5)[master_index]
prop.restore[,5] <- replace(v, which(is.na(v)), 0)
summary(prop.restore)

for (i in 1:dim(prop.restore)[1]){
  prop.restore[i,] <- prop.restore[i,] / sum(prop.restore[i,])
}
summary(prop.restore)
summary(rowSums(prop.restore))


# need to fill in NAs here based on ecoregion data
# r.ernatveg <- rcompose("./Ecoregion2017_natveg_classification/Ecoregions2017_Eco_valueNames_4.9km_Molweide_reclassificado_ESA.tif")
r.ernatveg <- rcompose("./5km/Ecoregions/Ecoregions2017__4.9km_Molweide_reclassificado.tif")

recs <- which(is.na(rowSums(prop.restore)))
length(recs)


# ignore this: just preparing a table for Carlos to work with:
# problem.coords <- xyFromCell(r.er, master_index[recs])
# png(file="problem_coordinates.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.er, axes=F)
# points(problem.coords, pch=19, col="black", cex=1)
# dev.off()
# write.csv(data.frame(masterindex=recs, problem.coords), file="problem_coords.csv", row.names=FALSE)

# scp 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/problem_c*.png .



# v.ernv <- values(r.ernatveg)[master_index[recs]]
# table(v.ernv)
# that provides no information for filling them in, so just allocate desert?

# these values need to be manually filled in with the patch data Carlos sends. 
# but for now I just make them all desert so we can at least run the code:
prop.restore[recs, 1] <- 0
prop.restore[recs, 2] <- 0
prop.restore[recs, 3] <- 1
prop.restore[recs, 4] <- 0
prop.restore[recs, 5] <- 0
summary(prop.restore)

save(prop.restore, file=paste0(dir, "prop.restore.RData"))






### SPEED UP EXTINCTION RISK CALCULATIONS WITH THESE DATA STRUCTURES

usphab_proc <- unique(sphabm_proc)
usphab_index <- list()
for (i in 1:dim(usphab_proc)[1]){
  ids <- c()
  for (j in 1:dim(sphabm_proc)[1]){
    if (identical(usphab_proc[i,], sphabm_proc[j,])) ids <- c(ids, j)
  }
  usphab_index[[i]] <- ids
  message(paste("length: ", length(ids)))
}

save(usphab_proc, file=paste0(dir, "usphab_proc.RData"))
save(usphab_index, file=paste0(dir, "usphab_index.RData"))



### SIMULATE COUNTRY LEVEL TARGETS

# how much agricultural and pasture land in each country
v.c <- values(r.cntry)[terrestrial_index]
v.ag <- values(r.plc1)[terrestrial_index]
v.pa <- values(r.plc2)[terrestrial_index]

cntry.area.ag <- rep(0, 241)
cntry.area.past <- rep(0, 241)

for (i in 1:241){
  recs <- which(v.c == i)
  cntry.area.ag[i] <- sum(v.ag[recs], na.rm=T)
  cntry.area.past[i] <- sum(v.pa[recs], na.rm=T)
}
summary(cntry.area.ag)
summary(cntry.area.past)


# get list of country ids that have cells in master_list
v <- values(r.cntry)[master_index]
table(v)


























# ### OTHER CODE FOR REPORTING / ANALYSIS ######################################################
# 
# # no need to pay attention to any of this
# 
# ### TOTAL AREA CALCULATIONS
# 
# # cell area in km^2
# A <- (res(r.plc1)[1])^2 / 1E6
# 
# totalareas <- rep(0, 10)
# totalareas[1] <- sum(values(r.plc1), na.rm=TRUE)
# totalareas[2] <- sum(values(r.plc2), na.rm=TRUE)
# totalareas[3] <- sum(values(r.plc3), na.rm=TRUE)
# totalareas[4] <- sum(values(r.plc4), na.rm=TRUE)
# totalareas[5] <- sum(values(r.plc5), na.rm=TRUE)
# totalareas[6] <- sum(values(r.plc6), na.rm=TRUE)
# totalareas[7] <- sum(values(r.plc7), na.rm=TRUE)
# totalareas[8] <- sum(values(r.plc8), na.rm=TRUE)
# totalareas[9] <- sum(values(r.plc9), na.rm=TRUE)
# totalareas[10] <- sum(values(r.plc10), na.rm=TRUE)
# totalareas <- totalareas * A / 1E6
# 
# # land area earth: 510.1 million km2
# sum(totalareas)
# 
# labels <- c("cropland", "cultivated grassland", "herbaceous cover", "forest", "wetlands", "desert", "natural_grasslands", "shrubland", "non-restorable", "water")
# 
# for (i in 1:length(labels)){
# 	cat(labels[i], "\t", format(totalareas[i], digits=3), "\n")
# }
# 
# 
# 
# ### NOTES ####################################################################
# 
# 
# # land cover codes
# 
# # 1 forest
# # 2 wetland
# # 3 desert
# # 4 crop
# # 5 natural grassland
# # 6 cultivated grassland
# # 7 shrubland
# # 8 not restorable
# # 9 water
# # 10 class 11 - herbaciouos cover
# 
# 
# # dummy variables for testing:
# 
# test.sphabm <- matrix(0, nrow=10, ncol=5)
# test.sphabm[1,1] <- 1
# test.sphabm[2,2] <- 1
# test.sphabm[3,3] <- 1
# test.sphabm[4,4] <- 1
# test.sphabm[5,5] <- 1
# test.sphabm[6,1] <- 1
# test.sphabm[6,2] <- 1
# test.sphabm[7,3] <- 1
# test.sphabm[7,5] <- 1
# test.sphabm[8,1] <- 1
# test.sphabm[8,2] <- 1
# test.sphabm[9,2] <- 1
# test.sphabm[9,3] <- 1
# test.sphabm[9,4] <- 1
# test.sphabm[10,] <- 1
# 
# test.sphabm
# 
# 
# 
# 
# # testing extinction risk code
# x <- seq(0, 1000, len=11)
# y <- extinction.risk(x, 1000)
# plot(x, y, type="l")
# 
# for (i in 1:length(x)){
# 	message(calc.extinction.slope(x[i], 1000))
# }
# 
# 
# 
# 
# # sum the species ranges
# v <- values(rcompose(spp_raster_names[1]))
# for (i in 2:length(spp_raster_names)){
# 	v <- v + values(rcompose(spp_raster_names[i]))
# }
# r.allmammals <- rcompose(spp_raster_names[1])
# values(r.allmammals) <- v
# 
# save(r.allmammals, file="r.allmammal.RData")
# 
# 
# # rename 's/^b_//' *.tif
# # rename -n 's/_moll//' *.tif
# 
# #rename 's/jpg?.*/jpg/' *
# 
# load(file=paste0(dir, "master_index.RData"))
# load(file=paste0(dir, "terrestrial_index.RData"))
# load(file=paste0(dir, "species_index_list_proc.RData"))
# load(file=paste0(dir, "species_index_list_terr_proc.RData"))
# load(file=paste0(dir, "spp_taxon_proc.RData"))
# 
# 
# v.m <- rep(0, length(master_index))
# v.a <- rep(0, length(master_index))
# for (i in 1:length(species_index_list_proc)){
# 	if (spp_taxon_proc[i] == "A"){
# 		v.a[species_index_list_proc[[i]]] <- v.a[species_index_list_proc[[i]]] + 1
# 	} else {
# 		v.m[species_index_list_proc[[i]]] <- v.m[species_index_list_proc[[i]]] + 1
# 	}
# }
# 
# r.amph <- r.terr
# values(r.amph) <- NA
# values(r.amph)[master_index] <- v.a
# save(r.amph, file="r.amph.RData")
# r.mamm <- r.terr
# values(r.mamm) <- NA
# values(r.mamm)[master_index] <- v.a
# save(r.mamm, file="r.mamm.RData")
# 
# 
# plot.pu.map(v.m, fname="map_mammals.png", bias=1)
# plot.pu.map(v.a, fname="map_amphibians.png", bias=1)
# 
# plot.pu.map(prop.crop + prop.cultg, fname="map_available_restoration.png", bias=1)
# 
# 
# 
# 
# # create map of cells to fill in opportunity cost layer
# load(file="terrestrial_lands.RData")
# 
# 
# v <- values(r.oc)[master_index]
# recs <- which(is.na(v))
# length(recs)
# 
# r.tmp <- r.terr
# values(r.tmp)[master_index[recs]] <- 2
# 
# png(file="oppcost_NA.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.tmp, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey70", "red"), axes=F, box=F, legend=F)
# dev.off()
# 
# 
# v <- values(r.ocg)[master_index]
# recs <- which(is.na(v))
# length(recs)
# 
# r.tmp <- r.terr
# values(r.tmp)[master_index[recs]] <- 2
# 
# png(file="oppcost_grass_NA.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.tmp, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey70", "red"), axes=F, box=F, legend=F)
# dev.off()
# 
# 
# v <- values(r.occ)[master_index]
# recs <- which(is.na(v))
# length(recs)
# 
# r.tmp <- r.terr
# values(r.tmp)[master_index[recs]] <- 2
# 
# png(file="oppcost_crop_NA.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.tmp, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey70", "red"), axes=F, box=F, legend=F)
# dev.off()
# 
# 
# v <- values(r.ocg)[master_index]
# recs <- which(v == 0)
# length(recs)
# 
# r.tmp <- r.terr
# values(r.tmp)[master_index[recs]] <- 2
# 
# png(file="oppcost_grass_zero.png", width=2000, height=1000)
# par(mar=c(0,0,0,0))
# plot(r.tmp, xlim=c(-1.2E7, 1.6E7), ylim=c(-6.5E6, 8828767), col=c("grey90", "grey70", "red"), axes=F, box=F, legend=F)
# dev.off()
# 
# 
# # plot of sum of extinction risks
# 
# 
# v <- rep(0, length(master_index))
# l <- rep(0, length(species_index_list_proc))
# 
# for (i in 1:length(species_index_list_proc)){
# 	l[i] <-  (length(species_index_list_proc[[i]]))
# 	#if (max(species_index_list_proc[[i]]) > length(master_index)) message(i)
# 	if (length(species_index_list_proc[[i]]) > 0){
# 		v[species_index_list_proc[[i]]] <- v[species_index_list_proc[[i]]] + exrisk.t0[i]
# 		if (length(v) > length(master_index)) {message(i); break}
# 		# if (max(species_index_list_proc[[i]]) > length(master_index)) message(i)
# 	} 
# }
# 
# 
# v <- rep(0, length(master_index)]
# for (i in 1:length(species_index_list_proc)){
# 	message(i)
# 	message(length(species_index_list_proc))
# 	if (max(species_index_list_proc[[i]]) > length(master_index)) message(i)
# 	v[species_index_list_proc[[i]]] <- v[species_index_list_proc[[i]]] + exrisk.t0[i]
# }
# 
# plot.pu.map(v, fname="biodiversity_sum_extinction_risk_map.png", bias=1)
# 
# 
# 
# 
# 
# 
# # check for duplicate IDs among taxons in the species table
# tab <- table(spid)
# recs <- which(tab > 1)
# tab[recs]
# dup <- as.numeric(names(tab)[recs])
# recs <- which(spid %in% dup)
# 
# recs <- which(spid == dup[1])
# spp_taxon[recs]
# 
# mix <- rep("", length(dup))
# for (i in 1:length(dup)){
# 	recs <- which(spid == dup[i])
# 	mix[i] <- paste(spp_taxon[recs], collapse="")
# }
# table(mix)
# 
# 
# # are the species IDs unique among taxons?
# uid <- unique(sphabdf$taxonid)
# for (i in 1:length(uid)){
# 	if (length(unique(sphabdf$class_name[which(sphabdf$taxonid == uid[i])])) > 1){
# 		message(paste("Not unique among taxons: ", uid[i]))
# 	}
# }
# 
# 
# 
# 
# # do the OA rcomposes sum to 1? (they should not)
# 
# v1 <- values(r.oa1)[master_index]
# v2 <- values(r.oa2)[master_index]
# v3 <- values(r.oa3)[master_index]
# v4 <- values(r.oa4)[master_index]
# v5 <- values(r.oa5)[master_index]
# vsum <- v1 + v2 + v3 + v4 + v5
# 
# # what about the 600+ species that have zero pre-settlement habitat?
# 
# load(file=paste0(dir, "species_cell_cnt.RData"))
# recs <- which(is.na(exrisk.t0))
# species_cell_cnt[recs]
# 
# table(sphabm[recs,1])
# table(sphabm[recs,2])
# table(sphabm[recs,3])
# table(sphabm[recs,4])
# table(sphabm[recs,5])
# 
# spid_proc[recs]
# table(spp_taxon_proc[recs])
# 
# eliminate_spid_0_habitat <- spid_proc[recs]
# save(eliminate_spid_0_habitat, file=paste0(dir, "eliminate_spid_0_habitat.RData"))
# 
# 
# 
# # scp 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/inputdata_v1/species_index_list.RData .
# 
# # scp countries-code.tif 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
# # scp opportunity_cost_4.9km.tif 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
# # scp ESA_landuse*.tif 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
# 
# # scp countries-code.tif 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
# # scp countries-code.tif 10.12.193.171:/home/hbeyer/Documents/proj/iis_global_spp/
# 
# 



