
# Code to produce country level summaries of benefits and costs for a solution.
# Also produces error estimates for those values. For carbon, this is based on
# the Delta C raster +/- 1.96 * SD, where SD is the standard deviation extracted
# from Cata's raster. Note that *negative values are permitted* in the Delta C
# raster, therefore it is acceptable for the error estimates to also include
# negative values.


library(raster)

# set the global input directory variable
setwd('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_priorization/')
dir <- "inputdata_v6/"

# load the master index:
load(paste0(dir, "master_index.RData"))



# constant to convert tonnes/ha to mega tonnes (also = peta grammes)
const_cb <- 1E-4

# constant to convert oc to $millions
const_oc <- 1E-6

# constant to convert area to millions of km2
const_area <- 1 #1E-6

# constant to convert mean abs extinction risk delta proportion to a percent
const_bd <- 1E2



### PREPROCESSING - ONLY NEED TO RUN THIS ONCE EVERY TIME THE MASTER INDEX CHANGES ###
# Creates a data structure that is saved to the input data folder and can be reused.

# load the CB and CB SD rasters
r.cb = raster('./rawdata/5km/others/DELTA_C_BiomassSoil30cm_v12.1.tif')
r.cbsd = raster('./rawdata/5km/others/C_SD_Biomass_Combined_COMPLETE_v12.1.tif')

# Check that projections and dimensions are identical:
r.cb
r.cbsd

# extract the CB SD values:
cbsd <- values(r.cbsd)[master_index]
summary(cbsd)
# fill in NA values with mean carbon
# this is not needed as there are currently no NA cells, but it is safe to run this code
# because it does nothing if there are no NA cells
if (length(which(is.na(cbsd))) > 0) cbsd <- replace(cbsd, which(is.na(cbsd)), mean(cbsd, na.rm=T))
save(cbsd, file=paste0(dir, "cbsd.RData"))
summary(cbsd)


# EXTRACT THE GRADIENT RASTER VALUES
r.gr <- raster('./display_results_WRLD_v8/scen_world-world_gradient-cb-bd-oc_total.restored.pu.tif')
r.gr

gradient <- values(r.gr)[master_index]
summary(gradient)
# fill in NA values with mean carbon
# this is not needed as there are currently no NA cells, but it is safe to run this code
# because it does nothing if there are no NA cells
if (length(which(is.na(gradient))) > 0) gradient <- replace(gradient, which(is.na(gradient)), mean(gradient, na.rm=T))
save(gradient, file=paste0(dir, "gradient.RData"))
summary(gradient)

### END PREPROCESSING ###




### FUNCTION THAT TAKES A SCENARIO SOLUTION AS INPUT AND PRODUCES COUNTRY-LEVEL SUMMARIES

# assumes that the global "dir" variable has been set that defines the input folder

# the function writes a csv file that can be opened in Excel in order to easily create
# a table that can be embedded in Word - so the countries are ordered alphabetically


# Alvaro: I set this as a variable because it is likely to be different on your system:
countrycodesfile = "../world-prod-estimates/countries-shp/countries-code.csv"


# LOAD THE FUNCTION BELOW NOW, THEN RUN IT:
#resultfile = "country-analysis/scen_cb-bd-oc_res.total.restored.pu_w_5.RData"
#outfile = "scen_cb-bd-oc_res.total.restored.pu_w_5.csv"
resultfile = './opt_results_CBD_v8/scen_cb-bd-oc_res.total.restored.pu_w_3.RData'
outfile = "./display_results_CBD_v8/scen_cb-bd-oc_res.total.restored.pu_w_3-per_coutry_post_processed.csv"

postproc.country.summary(resultfile, countrycodesfile, outfile)





postproc.country.summary <- function(resultfile, countrycodesfile, outfile){


	# load standard data objects from input data folder
	load(paste0(dir, "cntry.RData"))
	load(paste0(dir, "cb.RData"))
	load(paste0(dir, "cbsd.RData"))
	load(paste0(dir, "gradient.RData"))
	load(paste0(dir, "A.RData"))
	load(file=paste0(dir, "exrisk.t0.RData"))
    load(file=paste0(dir, "ocg.RData"))
    load(file=paste0(dir, "occ.RData"))
    load(file=paste0(dir, "prop.crop.RData"))
    load(file=paste0(dir, "prop.cultg.RData"))
    oc <- (prop.crop / (prop.crop + prop.cultg)) * occ + (prop.cultg / (prop.crop + prop.cultg)) * ocg 
    ub <- (prop.crop + prop.cultg)
	load(file=paste0(dir, "species_index_list_proc.RData"))
	load(file=paste0(dir, 'habarea.t0.RData'))
	load(file=paste0(dir, 'habarea.max.RData'))
	load(file=paste0(dir, "prop.restore.RData"))
	load(file=paste0(dir, "sphabm_proc.RData"))
	source('functions.r')
	cntrycodes <- read.csv(countrycodesfile)
	c.ord <- order(cntrycodes$NAME)

	load(resultfile)

	# calculate coutry level CB returns with upper and lower CI estimates
	nc <- length(c.ord)
	cname <- rep("", nc)
	sarea <- rep(NA, nc)
	savail <- rep(NA, nc)
	scb <- rep(NA, nc)
	scbll <- rep(NA, nc)
	scbul <- rep(NA, nc)
	sert0 <- rep(NA, nc)
	sert0n <- rep(NA, nc)
	serrest <- rep(NA, nc)
	exavoid <- rep(NA, nc)
	exavoidsd <- rep(NA, nc)
	mgr <- matrix(0, nrow=nc, ncol=10)
	soc <- rep(NA, nc)

	# priorty classes for gradient reporting:
	prcl <- c(0, 5, 10, 15, 20, 30, 40, 55, 70, 85, 100)

	for (i in 1:length(c.ord)){
		message(paste0("Country ", i, " of ", length(c.ord)))

		cname[i] <- cntrycodes$NAME[c.ord[i]]
		recs <- which(cntry == cntrycodes$CODE[c.ord[i]])
		if (length(recs) > 0){
			# CB
			sarea[i] <- sum(res.total.restored.pu[recs] * A * const_area)
			savail[i] <- sum(ub[recs] * A * const_area)
			scb[i] <- sum(cb[recs] * res.total.restored.pu[recs] * A * const_cb)
			scbll[i] <- sum((cb[recs] - 1.96 * cbsd[recs]) * res.total.restored.pu[recs] * A * const_cb)
			scbul[i] <- sum((cb[recs] + 1.96 * cbsd[recs]) * res.total.restored.pu[recs] * A * const_cb)

			# OC
			soc[i] <- sum(oc[recs] * res.total.restored.pu[recs] * A * const_oc)

			# BD
			# find all species with original ranges that overlap the country
			spp_idx <- c()
			for (j in 1:length(species_index_list_proc)){
				if (length(which(recs %in% species_index_list_proc[[j]])) > 0) spp_idx <- c(spp_idx, j)
			}

			# for each species, calculate extinction risk reduction resulting from in-country restoration
			sert0n[i] <- length(spp_idx)
			sert0[i] <- sum(exrisk.t0[spp_idx]) * 100 / length(spp_idx)

			# calculate reduction in extinction risk following in-country restoration

			# create a temporary vector representing restoration only within this country:
			cntry.restored <- rep(0, length(res.total.restored.pu))
			cntry.restored[recs] <- res.total.restored.pu[recs]

			# allocate that restoration to habitat types
			delta.vegtype.pu <- prop.restore * cntry.restored

			ns <- length(species_index_list_proc)

			# convert the change in habitat per pu to change in habitat per species
			delta.hab.spp <- rep(0, ns)
			for (k in 1:ns){
				for(j in 1:5){
					if (sphabm_proc[k,j] == 1){
						delta.hab.spp[k] <- delta.hab.spp[k] + sum(delta.vegtype.pu[species_index_list_proc[[k]], j], na.rm=T)
					}
				}
			}
			delta.hab.spp <- delta.hab.spp * A
			cntry.exrisk <- extinction.risk(habarea.t0 + delta.hab.spp, habarea.max, z=0.25)

			serrest[i] <- sum(cntry.exrisk[spp_idx]) * 100 / length(spp_idx)
			exavoid[i] <- sum(exrisk.t0[spp_idx] - cntry.exrisk[spp_idx])
			exavoidsd[i] <- sqrt(sum(exrisk.t0[spp_idx] * (1 - exrisk.t0[spp_idx])) + 
				sum(cntry.exrisk[spp_idx] * (1 - cntry.exrisk[spp_idx])))

			if (sarea[i] > 0){
				for (k in 1:10){
					recs2 <- which(gradient[recs] > prcl[k] & gradient[recs] <= prcl[k+1])
					if (length(recs2) > 0) mgr[i,k] <- (sum(res.total.restored.pu[recs[recs2]] * A * const_area) * 100) / sarea[i]
				}
			}

		}

	}


	df <- data.frame(cname, savail, sarea, scb, scbll, scbul, sert0n, sert0, serrest, exavoid, exavoidsd, soc, mgr)
	names(df) <- c("Country", "Area available (km2)", "Area restored (km2)", "Delta C (Mt)", 
		"Delta C lower (Mt)", "Delta C upper (Mt)", "N species", "Mean ext. risk (pre-restoration; %)",
		 "Mean ext. risk (post-restoration; %)", "Extinctions avoided (mean)", 
		 "Extinctions avoided (SD)", "Opportunity cost ($M)", "Gradient scenario (%; 0-5%)",
		 "Gradient scenario (%; 5-10%)", "Gradient scenario (%; 10-15%)", "Gradient scenario (%; 15-20%)",
		 "Gradient scenario (%; 20-30%)", "Gradient scenario (%; 30-40%)", "Gradient scenario (%; 40-55%)",
		 "Gradient scenario (%; 55-70%)", "Gradient scenario (%; 70-85%)", "Gradient scenario (%; 85-100%)")

	# remove records with NA:
	df <- df[-which(is.na(sarea)),]

	write.csv(df, file=outfile)
	return(1)
}









