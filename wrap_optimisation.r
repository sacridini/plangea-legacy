# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

# This is the starting point of the analyses in this version.

# It divides the processing in blocks of scenario runs, defined by the amount of
# restoration allowed in each PU. Each one of those blocks have control
# parameters defining the range of scenarios to be included (more below), but in
# general each block represents a full run of the previous version, in a given
# setting of overall limit to the restoration of each PU.

# Having defined the control parameters for a given scenario block, the code
# goes on to call optimisation.r, which was the starting point of the previous
# version. In total it calls optimisation.r four times -- three for different
# restoration limits, and a last time for building a gradient map of priority
# areas for restoration across all the available area.

# New objects defined in this script:
# info (logical): determines if plots should be printed with info on their title
# print.steps (logical): should plot titles contain step information
# ub.perc.constraint: vector ub feeding into solver will be multiplied by this
# ublim.suffix (string): used in naming output files for that scenario block

# The following objects determine, for each scenario block, the range of
# scenarios parameters to be included in the run of the block. The scenario
# parameters are defined in wrap_scenarios.r, and all of them span one of the
# three ranges (which will only be properly defined in wrap_scenario.r):

# scc.range: 1: no country constraints sent to solver
#            2: with country constraints
# sct.range: 1: Bonn-challenge restoration target (1.5 million sq.km)
#            2: NY declaration restoration target (3.5 million sq.km)
#            3: sum of sparable land in each country minus demanded land
#            4: CBD target (15% of restorable area)
#            5: double of CBD target
#            6: all the restorable area (useful for gradient maps)
# scb.range: 1: "scen_cb-oc" benchmark scenario included
#            2: "scen_cb" benchmark scenario included
#            3: "scen_bd-oc"  benchmark scenario included
#            4: "scen_bd" benchmark scenario included
#            5: "scen_cb-bd-oc" benchmark scenario included
#            6: "scen_cb-bd" benchmark scenario included
#            7: "scen_oc" benchmark scenario included
#            8: "scen_rnd" benchmark scenario included
#            9: "scen_world" ("cb-bd-oc" with steps to create gradient map)

# The control objects defined in each scenario block work as follows:
# ublim.cty.range: controls which elements in scc.range to be included
# target.range: controls which elements in sct.range to be included
# bench.range: controls which elements in scb.range to be included

# Beware that large targets in sct.range might not be achievable in scenario
# runs with low values of ub.perc.constraint

# Also, the changes I needed to do in order to integrate the previous version
# of the code to all the scenario runs, especially in optimisation.r and
# run_optimisation.r was rendering wrong results on the random allocation
# scenarios ("scen_rnd"). I made many corrections to the code since the last
# time I tried to run "scen_rnd", but I haven't checked if they corrected the
# random allocation scenario as well.


# PRE-PROCESSING (INCLUDING AUXILIARY ANALYSES: OA AND OPPORTUNITY COSTS) ######
#source('preprocessing.r', echo=T)
#rm(list=ls(all=T))
#Sys.sleep(1)
#gc()


# SCENARIO BLOCK 1: OVERALL-RESTORATION LIMITS ONLY ############################
info = F # plangea: not needed
print.steps = F # plangea: not needed
z.sd = 0.1 # plangea: config.json -> $variables$calc_bd$z_std_dev
DR.sd = 0.05 # plangea: too specific of how OC was implemented on legacy
PrC.relative.var = 0.25 # plangea: too specific of how OC was implemented on legacy
PrG.relative.var = 0.25 # plangea: too specific of how OC was implemented on legacy
quad.sd = T # plangea: quadratic propagation should be default
CL.prt = T # plangea: not needed
ub.perc.constraint = 1
ub.perc.const = c(0.25, 0.5, 0.75) # plangea: config.json -> $scenarios$upper_bound_limits
econ.ub = F # plangea: analysis using econometric limits to ub, specific of legacy
econ.ctrylim = F # plangea: analysis using econometric limits to ctry.lims, specific of legacy
ublim.suffix = ifelse(econ.ub, '-econ-ub_',  '') 
ublim.cty.range = 1 # plangea: config.json -> $scenarios$include_subregion_scenarios
target.range = 4 # plangea: outdated control mechanism of legacy
bench.range = 1:7 # plangea: outdated control mechanism of legacy
overwrite.nsteps = 5 # plangea: outdated control mechanism of legacy
wgt.range = 1 # plangea: outdated control mechanism of legacy
source("optimisation.r")

 for (ub.perc.constraint in seq_along(ub.perc.const)) {
   ublim.suffix = paste0('-econ-ub_', ub.perc.constraint)
   source("optimisation.r")
}

rm(list=ls(all=T))
Sys.sleep(1)
gc()



# SCENARIO BLOCK 2: COUNTRY LIMITS TO RESTORATION ##############################
info = T
print.steps = T
z.sd = 0.1
DR.sd = 0.05
PrC.relative.var = 0.25
PrG.relative.var = 0.25
quad.sd = T
CL.prt = T
wrld.form = 1
ub.perc.constraint = 1
econ.ub = F
econ.ctrylim = T
#flat.ctrylim.vals = seq(from=0.15, to=0.45, by=0.1)
flat.ctrylim.vals = 1 #plangea: config.json -> $scenarios$sub-region_scenarios$sub-region_flat_targets_multiplier
ublim.cty.range = 2
target.range = 4
bench.range = 1:7
overwrite.nsteps = 5
wgt.range = 1

for (flat.ctrylim in flat.ctrylim.vals){
  econ.ctrylim.suffix = ifelse(econ.ctrylim, '-econ-ctrylims_' ,'')
  ublim.suffix = ifelse(exists('flat.ctrylim'),
                        paste0('-ctrylim_', flat.ctrylim, econ.ctrylim.suffix),
                        econ.ctrylim.suffix)
  ublim.suffix = paste0(ifelse(econ.ub, '-econ-ub' , ''),ublim.suffix)
  source("optimisation.r")
}
#  if (exists('flat.ctrylim')){rm(flat.ctrylim)}
  
#  flat.ctrylim.df = c()
#  for (ns in 1:nsteps){
    #  load(paste0(outdir, scen, "_step.res_", ns, ".RData"))
    #  step.ras = step.ras + ((nsteps-i) * step.res)
#    flat.ctrylim.df = rbind(flat.ctrylim.df,
#                            postprocess.grad(dir, outdir, ns=ns,
#                                             filename=paste0(scen, "_res.total.restored.pu_step_", ns, ".RData"#)))
#  }
#  write.csv(flat.ctrylim.df, file=paste0("./display_results_CBD_v8/world_gradient_results",ublim.suffix, ".csv"), row.names=F)
#}

rm(list=ls(all=T))
Sys.sleep(1)
gc()



# SCENARIO BLOCK 3: STEP-WISE GLOBAL RESTORATION GRADIENT ######################
info = T
print.steps = T
z.sd = 0.05
DR.sd = 0.01
PrC.relative.var = 0.25
PrG.relative.var = 0.25
quad.sd = T
CL.prt = T
ub.perc.constraint = 1
econ.ub = F
econ.ctrylim = F
wrld.suffix = c('cb-bd') # , 'cb-bd-oc', 'bd', 'bd-oc', 'cb', 'cb-oc', 'oc')
wgt.range = 5
wrld.res.df = c()

for (wrld.loop in wrld.suffix){
  ublim.suffix = paste0('-world_gradient-', wrld.loop)
  ublim.suffix = paste0(ifelse(econ.ub, '-econ-ub' , ''),ublim.suffix)
  wrld.form = which(wrld.suffix == wrld.loop)
  ublim.cty.range = 1
  target.range = 6
  bench.range = 9
  source("optimisation.r")
  
  #load(paste0(outdir, scen, "_step.res_1.RData"))
  #step.ras = step.res
  
  for (wgt in wgt.range){
    for (ns in 1:nsteps){
      #  load(paste0(outdir, scen, "_step.res_", ns, ".RData"))
      #  step.ras = step.ras + ((nsteps-i) * step.res)
      wrld.res.df = rbind(wrld.res.df,
                          postprocess.grad(dir, outdir, ns=ns, filename=paste0(scen, '_res.total.restored.pu_w', wgt,'_step_', ns, '.RData')))
    }    
  }

  
  #plot.pu.map(step.ras[master_index], fname=paste0('./opt_results_WRLD_v8/scen_world_incremental-', ublim.suffix,'.png'))
  
}

write.csv(wrld.res.df, file="./display_results_WRLD_v8/world_gradient_results.csv", row.names=F)

rm(list=ls(all=T))
Sys.sleep(1)
gc()



# SCENARIO BLOCK 4: UPPER BOUNDS LIMITED TO FRACTION OF THE CELL ###############
info = F # Should plots be printed with diagnostics info?
print.steps = F
z.sd = 0.1
DR.sd = 0.05
PrC.relative.var = 0.25
PrG.relative.var = 0.25
quad.sd = T
CL.prt = T
ub.perc.vals = seq(from=0.15, to=0.95, by=0.1)
econ.ub = F
econ.ctrylim = F
ublim.cty.range = 1
target.range = 4
bench.range = 1:7
overwrite.nsteps = 1
#wgt.range = 1

for (ub.perc.constraint in ub.perc.vals){
    ublim.suffix = paste0('-ublim_',round(ub.perc.constraint,2))
    ublim.suffix = paste0(ifelse(econ.ub, '-econ-ub' , ''),ublim.suffix)
    source("optimisation.r")
}

rm(list=ls(all=T))
Sys.sleep(1)
gc()



# SCENARIO BLOCK 5: GLOBAL RESTORATION TARGET MATCHING UB LIMIT ################
info = F # Should plots be printed with diagnostics info?
print.steps = F
z.sd = 0.1
DR.sd = 0.05
PrC.relative.var = 0.25
PrG.relative.var = 0.25
quad.sd = T
CL.prt = T
# ub.perc.vals = seq(from=0.05, to=0.15, by=0.05)
ub.perc.vals = seq(from=0.15, to=1, by=0.1)
econ.ub = F
econ.ctrylim = F
ublim.cty.range = 1
target.range = 4
wgt.range = c(1)


match.flag = F # defines if ub limit should be matched to global limit (T) or not (F)

if (match.flag){ # Matching ub limit to global (tgt.mod) limit
  bench.range = 5
  overwrite.nsteps = 1
  for (ub.perc.constraint in ub.perc.vals){
    tgt.mod.generalized = ub.perc.constraint
    ublim.suffix = paste0('-ublim-matching-global-at_',round(ub.perc.constraint,2))
    source("optimisation.r")
    if(!exists('iter.results')){iter.results = results.df}
    if(as.character(iter.results[1,1]) != as.character(results.df[1,1])){iter.results = rbind(iter.results, results.df)}
  }
  write.csv(iter.results, file="./display_results_CBD_v8/ublim-matching-global_combined-results.csv", row.names=F)  
} else {# Changing global (tgt.mod) limit keeping ub limit equal to 1
  for (tgt.mod.generalized in ub.perc.vals){
    bench.range = 5 #c(2,4,5,7)
    overwrite.nsteps = 5
    ublim.suffix = paste0('-global-limit-at_',round(tgt.mod.generalized,2))
    ub.perc.constraint = 1
    source("optimisation.r")
    if(!exists('iter.results')){iter.results = results.df}
    if(as.character(iter.results[1,1]) != as.character(results.df[1,1])){iter.results = rbind(iter.results, results.df)}
  }
  #write.csv(iter.results, file="./display_results_CBD_v8/global_limits-combined-results.csv", row.names=F)
  #write.csv(iter.results, file="./display_results_CBD_v8/compromise_tradeoff-combined-results.csv", row.names=F)
}
