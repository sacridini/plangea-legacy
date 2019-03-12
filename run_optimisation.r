# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

# This is the general purpose code that gets run for all scenarios.
# It is important not to call this from within a looping structure that uses indices
# implemented in this script: s, i, j

# The only addition to this version was a code block defining the optimisation
# problem, and doing a call to Rsymphony instead of gurobi. The block is
# self-contained, and clearly defined right after the original call to gurobi.
# To resume the call to gurobi one must simply de-comment the original call
# and comment the RSymphony block.

# results data objects:
res.prop.restored.pu <- matrix(0, nrow=np, ncol=nsteps)
res.total.restored.pu <- rep(0, np)
res.area.restored.spp <- matrix(0, nrow=ns, ncol=nsteps)
res.total.restored.spp <- rep(0, ns)
res.exrisk <- matrix(0, nrow=ns, ncol=nsteps)
res.exrisk.sd <- matrix(0, nrow=ns, ncol=nsteps)
res.objval <- rep(0, nsteps)

res.grad = res.total.restored.pu

delta.hab.pu <- rep(0, np) #<- matrix(0, nrow=np, ncol=5)


# loop over restoration increments
for (s in 1:nsteps){

  print(paste0('Running scenario ', save.name, ' with target ', sct,
               ' (', round(restoration.area, digits=2), ' sq.km). Step: ', s,
               '/', nsteps, '. Weights ', w, ' (', w.cb, ' CB, ', w.bd,
               ' BD). Time elapsed: ', format(round(Sys.time() - start,1))))
  
	# set up Gurobi model in R
	model <- list()
	model$obj <- objfun
	model$modelsense <- "max"
	model$vtype <- "S"
	model$lb <- rep(0, np)
	model$ub <- ub
	model$A <- constr
	model$rhs <- rhs
	model$sense <- sense
	params <- list(Threads=4) #Presolve=2, MIPGap=0.005,
	#result <- gurobi(model,params)
	
	# R_SYMPHONY ALTERNATIVE BLOCK ###############################################
	model$ub[model$ub < model$lb] = model$lb[model$ub < model$lb]
	rsymphony_bounds = list(lower = list(ind=1:np, val=model$lb),
	                        upper = list(ind=1:np, val=model$ub))
	result = Rsymphony_solve_LP(obj = model$obj, mat = model$A,
	                            dir = model$sense, rhs = model$rhs,
	                            bounds = rsymphony_bounds, max=T)
	names(result)[1] = 'x'
	if (result$status == 0){result$status = "OPTIMAL"}
	##############################################################################

	if (result$status == "OPTIMAL"){

		res.objval[s] <- result$objval

		# record the restoration
		res.prop.restored.pu[,s] <- result$x
		res.total.restored.pu <- rowSums(res.prop.restored.pu) #res.total.restored.pu + result$x
		
		if (nsteps>1){
		  #res.grad[res.total.restored.pu == 0] = 1/nsteps
		  #res.grad[result$x > 0] = (2/nsteps) + (1 - 2/nsteps) * (1 + nsteps - s)
		  res.grad[result$x > 0] = 1 + nsteps - s
		  } else {res.grad = res.total.restored.pu}
		
		if (info) {
		  fname = NULL
		  if (print.steps){fname=paste0(outdir, scen, "_res.prop.restored.pu_step_", s, ".png")}
		  step.res = plot.pu.map(res.prop.restored.pu[,s], bias=1, fname = fname,
		                       info=T, area.print=T, step.count=T,
		                       localize = print.steps)
		  save(step.res, file=paste0(outdir, scen, "_step.res_", s, ".RData"))
		  save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_step_", s, ".RData"))
		  
		  if (exists('flat.ctrylim')){
		    ctry.lims = cbind(ctry.lims, sapply(world.csv$CODE, function(x){sum(model$A[x,] * (res.total.restored.pu) * A)}))
		    names(ctry.lims)[length(ctry.lims[1,])] = paste0('tot.res.A.s', s) }
		    #print(cbind(model$rhs[33]*nsteps, ctry.lims[33,-(1:2)]/A))
		  }

		# update the habitat proportions
		# be careful with the matrix multiplications here... easy to mess this up

		# allocate that restoration to habitat types
		delta.vegtype.pu <- prop.restore * res.total.restored.pu

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
		res.area.restored.spp[,s] <- delta.hab.spp
		res.total.restored.spp <- res.total.restored.spp + delta.hab.spp

		# recalculate bd weights based on revised habitat
		exrisk.slope <- calc.extinction.slope(habarea.t0 + delta.hab.spp, habarea.max, z=0.25)
		#exrisk.ts <- extinction.risk(habarea.t0 + delta.hab.spp, habarea.max, z=0.25)
		res.exrisk[,s] <- extinction.risk(habarea.t0 + delta.hab.spp, habarea.max, z=0.25)
		res.exrisk.sd[,s] <- extinction.risk.sd(habarea.t0 + delta.hab.spp, habarea.max, z=0.25, z.sd=z.sd)


		# reclaculate bd benefit
		bd.s <- calc.bd(exrisk.slope)
		#bd.s <- calc.bd(res.exrisk[,s])  # <------------------ testing code only

		# update the objective function
		objfun <- calc.objective.function(objfuncform, bd.s, cb, oc, w.cb, w.bd, wrld.form=wrld.form)

		# update the ub for each cell - this is the mechanism by which we restore incrementally
		# more habitat in each step
		ub <- ub - result$x

		if (g_diagnostic_plots){
			plot.pu.map(objfun, fname=paste0(outdir, scen, "_map.objfun_t", s, ".png"))
			plot.pu.map(res.total.restored.pu / (prop.crop + prop.cultg), fname=paste0(outdir, scen, "_map.total.restored.pu.pct_", w, "_t", s, ".png"))
		}
	} else {
		message("#### SOLUTION FAILURE ###########################")
		message(result$status)
		g_fatal_error <- TRUE
		break
	}

	# as a precaution, explicitly delete these variables to prevent inappropriate re-use
	if (exists("model")) rm(model)
	if (exists("result")) rm(result)
	if (exists("delta.hab.spp")) rm(delta.hab.spp)
	if (exists("delta.vegtype.pu")) rm(delta.vegtype.pu)
	if (exists("exrisk.slope")) rm(exrisk.slope)
	if (exists("bd.s")) rm(bd.s)

} # s

# Reset upper boundaries after step-s loop is finished
ub = (gap.agr * prop.crop + gap.grs * prop.cultg) * ub.perc.constraint

#save(res.objval, file=paste0(outdir, scen, "_res.objval_w_", w, ".RData"))
#save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
#save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
#save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
# save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
#save(res.exrisk, file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
