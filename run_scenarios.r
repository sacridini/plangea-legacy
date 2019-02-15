
### SCENARIOS #######################################################################

# note that for many of the scenarios w.cb and w.bd are not used, but we just set them to 1
# as a precautionary measure

### SCENARIO: CARBON ONLY #######################################################

# scenario name
scen <- "scen_cb-oc"

# constraint matrix
nconstr <- 1
constr <- matrix(0, nrow=nconstr, ncol=np)
sense <- rep("<", nconstr)
rhs <- rep(0, nconstr)

constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# multiple steps are only needed for BD
nsteps <- 1
rhs <- restoration.rhs

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

# calculate initial bd value - we reset this every time as a precautionary measure
#bd <- calc.bd(exrisk.t0.slope)

# set the objective function
objfuncform <- "cb-oc"
w.cb = 1
w.bd = 1
objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
# objfun <- cb / oc
# objective.function <- function(){ return(cb / oc) }

# there are no weights in this scenario but we still need to set w to something for file naming purposes
w <- "NA"

# plot map of objective function
plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_", w, ".png"))

# once all parameters are set, run the optimiser
source("run_optimisation.r")

# plot map of solution
plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
ub <- (prop.crop + prop.cultg)
plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))



### scenario name
scen <- "scen_cb"

# constraint matrix
nconstr <- 1
constr <- matrix(0, nrow=nconstr, ncol=np)
sense <- rep("<", nconstr)
rhs <- rep(0, nconstr)

constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# multiple steps are only needed for BD
nsteps <- 1
rhs <- restoration.rhs

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

# calculate initial bd value - we reset this every time as a precautionary measure
#bd <- calc.bd(exrisk.t0.slope)

# set the objective function
objfuncform <- "cb"
w.cb = 1
w.bd = 1
objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
# objfun <- cb
# objective.function <- function(){ return(cb) }

# there are no weights in this scenario but we still need to set w to something for file naming purposes
w <- "NA"

# plot map of objective function
plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_", w, ".png"))

# once all parameters are set, run the optimiser
source("run_optimisation.r")

# plot map of solution
plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
ub <- (prop.crop + prop.cultg)
plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))




### SCENARIO: BIODIVERSITY ONLY ###########################################

# scenario name
scen <- "scen_bd-oc"

# constraint matrix
nconstr <- 1
constr <- matrix(0, nrow=nconstr, ncol=np)
# sense <- rep("<=", nconstr)
# rhs <- rep(0, nconstr)

constr[1,] <- 1 * g_scalar_area
sense <- c("<=")


# multiple steps are needed for BD
nsteps <- 20
rhs <- restoration.rhs / nsteps

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

# calculate initial bd value
#bd <- calc.bd(exrisk.t0.slope)

# set the objective function
objfuncform <- "bd-oc"
w.cb = 1
w.bd = 1
objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
# objfun <- bd / oc
# objective.function <- function(){ return(bd / oc) }

# there are no weights in this scenario but we still need to set w to something for file naming purposes
w <- "NA"

# plot map of objective function
plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_t0.png"))

# once all parameters are set, run the optimiser
source("run_optimisation.r")

# plot map of solution
plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
ub <- (prop.crop + prop.cultg)
plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))



### scenario name
scen <- "scen_bd"

# constraint matrix
nconstr <- 1
constr <- matrix(0, nrow=nconstr, ncol=np)
# sense <- rep("<=", nconstr)
# rhs <- rep(0, nconstr)

constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# multiple steps are needed for BD
nsteps <- 20
rhs <- restoration.rhs / nsteps

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

# calculate initial bd value
#bd <- calc.bd(exrisk.t0.slope)

# set the objective function
objfuncform <- "bd"
w.cb = 1
w.bd = 1
objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
# objfun <- bd
# objective.function <- function(){ return(bd) }

# there are no weights in this scenario but we still need to set w to something for file naming purposes
w <- "NA"

# plot map of objective function
plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_t0.png"))

# once all parameters are set, run the optimiser
source("run_optimisation.r")

# plot map of solution
plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
ub <- (prop.crop + prop.cultg)
plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))





### SCENARIO: WEIGHTED TRADEOFF BETWEEN CARBON AND BIODIVERSITY ###################

# scenario name
scen <- "scen_cb-bd-oc"

# weights for carbon and biodiversity
# weightm <- matrix(1, nrow=9, ncol=2)
# weightm[1,] <- c(1,0)
# weightm[2,] <- c(1,0.1)
# weightm[3,] <- c(1,0.5)
# weightm[4,] <- c(1,1.5)
# weightm[5,] <- c(1,4)
# weightm[6,] <- c(1,10)
# weightm[7,] <- c(1,50)
# weightm[8,] <- c(1,100)
# weightm[9,] <- c(0,1)


# weights for carbon and biodiversity
weightm <- matrix(1, nrow=10, ncol=2)
weightm[1,] <- c(1,0)
weightm[2,] <- c(1,1)
weightm[3,] <- c(1,4)
weightm[4,] <- c(1,10)
weightm[5,] <- c(1,40)
weightm[6,] <- c(1,100)
weightm[7,] <- c(1,400)
weightm[8,] <- c(1,1000)
weightm[9,] <- c(1,4000)
weightm[10,] <- c(0,1)

save(weightm, file=paste0(outdir, scen, "_weightm.RData"))

########################################################################
# alternative weighting schemes for testing purposes
# # weights for carbon and biodiversity
# weightm <- matrix(1, nrow=3, ncol=2)
# weightm[1,] <- c(1,0)
# weightm[2,] <- c(1,100)
# weightm[3,] <- c(0,1)

# save(weightm, file=paste0(outdir, scen, "_weightm.RData"))

# # weights for carbon and biodiversity
# weightm <- matrix(1, nrow=6, ncol=2)
# weightm[1,] <- c(1,0)
# weightm[2,] <- c(1,0.0001)
# weightm[3,] <- c(1,0.001)
# weightm[4,] <- c(1,0.01)
# weightm[5,] <- c(1,0.1)
# weightm[6,] <- c(0,1)

# save(weightm, file=paste0(outdir, scen, "_weightm.RData"))
########################################################################


# constraint matrix
nconstr <- 1
# sense <- rep("<", nconstr)
# rhs <- rep(0, nconstr)
constr <- matrix(0, nrow=nconstr, ncol=np)
constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# multiple steps are needed for BD
nsteps <- 5
rhs <- restoration.rhs / nsteps


for (w in 1:dim(weightm)[1]){

	# essential to reset bd and ub before starting the next loop
	# calculate initial bd value
	#bd <- calc.bd(exrisk.t0.slope)

	# upper bound cannot exceed area available for restoration:
	ub <- (prop.crop + prop.cultg)

	# set the objective function
	objfuncform <- "cb-bd-oc"
	w.cb = weightm[w,1]
	w.bd = weightm[w,2]
	objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
	# objfun <- (weightm[w,1] * cb + weightm[w,2] * bd) / oc
	# objective.function <- function(){ return((weightm[w,1] * cb + weightm[w,2] * bd) / oc) }

	# plot map of objective function
	plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_", w, "_t0.png"))

	# once all parameters are set, run the optimiser
	source("run_optimisation.r")

	# plot map of solution
	plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
	ub <- (prop.crop + prop.cultg)
	plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))

}




### scenario name
scen <- "scen_cb-bd"

# weights for carbon and biodiversity
# weightm <- matrix(1, nrow=10, ncol=2)
# weightm[1,] <- c(1,0)
# weightm[2,] <- c(1,100)
# weightm[3,] <- c(1,500)
# weightm[4,] <- c(1,2500)
# weightm[5,] <- c(1,5000)
# weightm[6,] <- c(1,15000)
# weightm[7,] <- c(1,50000)
# weightm[8,] <- c(1,100000)
# weightm[9,] <- c(1,1000000)
# weightm[10,] <- c(0,1)
weightm <- matrix(1, nrow=9, ncol=2)
weightm[1,] <- c(1,0)
weightm[2,] <- c(1,0.1)
weightm[3,] <- c(1,0.5)
weightm[4,] <- c(1,1.5)
weightm[5,] <- c(1,4)
weightm[6,] <- c(1,10)
weightm[7,] <- c(1,50)
weightm[8,] <- c(1,100)
weightm[9,] <- c(0,1)

save(weightm, file=paste0(outdir, scen, "_weightm.RData"))


# constraint matrix
nconstr <- 1
# sense <- rep("<", nconstr)
# rhs <- rep(0, nconstr)
constr <- matrix(0, nrow=nconstr, ncol=np)
constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# here we are interested in a series of area targets so we can implement that by setting rhs to be the increment and adjusting nsteps
# targets should be translated to units of cells rather than area

# multiple steps are needed for BD
nsteps <- 5
rhs <- restoration.rhs / nsteps

for (w in 1:dim(weightm)[1]){

	# essential to reset bd and ub before starting the next loop
	# calculate initial bd value
	#bd <- calc.bd(exrisk.t0.slope)

	# upper bound cannot exceed area available for restoration:
	ub <- (prop.crop + prop.cultg)

	# set the objective function
	objfuncform <- "cb-bd"
	w.cb = weightm[w,1]
	w.bd = weightm[w,2]
	objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)

	# plot map of objective function
	plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_", w, "_t0.png"))

	# once all parameters are set, run the optimiser
	source("run_optimisation.r")

	# plot map of solution
	plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
	ub <- (prop.crop + prop.cultg)
	plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))

}




### SCENARIO: RANDOM ALLOCATION ####################################################

# note that here we do not need to use the scaled areas we use when calling Gurobi

# repeat this niters times
niters <- 10

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

for (k in 1:niters){

	scen <- paste0("scen_rnd", k)
	w <- "NA"
	nsteps <- 1

	recs <- random.allocation(restoration.area, ub * A)

	result <- list()
	x <- rep(0, length(ub))
	x[recs] <- ub[recs]
	result$x <- x

	# create the same data structures as are output from the optimisations
	res.prop.restored.pu <- matrix(0, nrow=np, ncol=nsteps)
	res.total.restored.pu <- rep(0, np)
	res.area.restored.spp <- matrix(0, nrow=ns, ncol=nsteps)
	res.total.restored.spp <- rep(0, ns)
	res.exrisk <- matrix(0, nrow=ns, ncol=nsteps)
	res.objval <- rep(0, nsteps)
	delta.hab.pu <- rep(0, np) 

	delta.hab.pu <- result$x
	# allocate that restoration to habitat types
	delta.vegtype.pu <- prop.restore * delta.hab.pu

	# convert the change in habitat per pu to change in habitat per species
	delta.hab.spp <- rep(0, ns)
	for (i in 1:ns){
		for(j in 1:5){
			if (sphabm_proc[i,j] == 1){
				delta.hab.spp[i] <- delta.hab.spp[i] + sum(delta.vegtype.pu[species_index_list_proc[[i]], j], na.rm=T)
			}
		}
	}
	delta.hab.spp <- delta.hab.spp * A
	res.area.restored.spp[,1] <- delta.hab.spp
	res.total.restored.spp <- res.total.restored.spp + delta.hab.spp

	res.exrisk[,1] <- extinction.risk(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25)
	res.prop.restored.pu[,1] <- result$x
	res.total.restored.pu <- result$x

	save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
	save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
	save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
	save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
	save(res.exrisk, file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))

	# plot map of solution
	plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
	ub <- (prop.crop + prop.cultg)
	plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))

}



### SCENARIO: OC ONLY ####################################################

# scenario name
scen <- "scen_oc"

# constraint matrix
nconstr <- 1
constr <- matrix(0, nrow=nconstr, ncol=np)
sense <- rep("<", nconstr)
rhs <- rep(0, nconstr)

constr[1,] <- 1 * g_scalar_area
sense <- c("<=")

# multiple steps are only needed for BD
nsteps <- 1
rhs <- restoration.rhs

# upper bound cannot exceed area available for restoration:
ub <- (prop.crop + prop.cultg)

# calculate initial bd value - we reset this every time as a precautionary measure
#bd <- calc.bd(exrisk.t0.slope)

# set the objective function
objfuncform <- "oc"
w.cb = 1
w.bd = 1
objfun <- calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd)
# objfun <- -oc + 130
# objective.function <- function(){ return(oc) }

# there are no weights in this scenario but we still need to set w to something for file naming purposes
w <- "NA"

# plot map of objective function
plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_", w, ".png"))

# once all parameters are set, run the optimiser
source("run_optimisation.r")

# plot map of solution
plot.pu.map(res.total.restored.pu, fname=paste0(outdir, scen, "_map.total.restored.pu_", w, ".png"))
ub <- (prop.crop + prop.cultg)
plot.pu.map(res.total.restored.pu / ub, fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, ".png"))

