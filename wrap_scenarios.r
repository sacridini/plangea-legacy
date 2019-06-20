# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

# Wrapper for calls to run_scenarios.r (itself a script connecting the new
# scenario parameter objects defined in wrap_optimisation.r and in here with the
# run_optimisation.r script from the original version, which actually gets to
# call the solver).

# In here the scenario parameters are defined, and each scenario is run from
# inside a loop spanning the three parameter ranges: target (sct.range), country
# constraints (scc.range), and benchmark scenarios (scb.range).

#source('run_scenario_fun.r')
#library('bit')

# Target scenarios (sct) #######################################################
tgt.mod = 0.15
if (exists('tgt.mod.generalized')){tgt.mod = tgt.mod.generalized}

rest.area = list(BONN=1.5e6, NYDC=3.5e6, CTRY=13.3719e6,
                 CBD = tgt.mod*sum((prop.crop + prop.cultg)*A),
                 CBD2 = 2*tgt.mod*sum((prop.crop + prop.cultg)*A),
                 WRLD = sum((prop.crop + prop.cultg)*A))

# Load large country.coefs matrix only if needed
country.coefs = c()
if (2 %in% ublim.cty.range){
  load(paste0(dir,'country.coefs.RData'))
  world.csv = read.csv(paste0(dir, 'countries-code.csv'))
  }

# Contraint equation coefficients
constr.list = list(unconstrained = matrix(g_scalar_area, nrow=1, ncol=np),
                   country.limits = country.coefs * g_scalar_area)

# Reading country-level limits to restoration
ctry.lims = read.csv(paste0(dir, 'restoration-constraints-per-country.csv'))

# Shaving-off unnecessary data from the table
ctry.lims = ctry.lims[-length(ctry.lims[,1]),c(1,length(ctry.lims[1,]))]

# Correcting country limits units from ha to sq.km
ctry.lims$total =  ctry.lims$total / 100 

if(!exists('ub.perc.constraint')){ub.perc.constraint=1}
if (econ.ctrylim){ctry.lims$total = econ.ctrylims$SQKM * ub.perc.constraint}

# Builds country limits depending on flat.ctrylim from wrap_optimisation
if (exists('flat.ctrylim')){ctry.lims$total = flat.ctrylim * ctry.lims$total}
  #ctry.lims$total = sapply(world.csv$CODE, function(x){flat.ctrylim * sum(country.coefs[x,] * (gap.agr*prop.crop + gap.grs*prop.cultg) * A)})
#} else {
  # Computing overall sparable land vs demanded land
#  dem.lnd = sum(ctry.lims$total[ctry.lims$total>0])
#  spa.lnd = sum(ctry.lims$total[ctry.lims$total<=0])
  
  # Correcting limits to restoration to account for compensation of demanded area
#  ctry.lims$total[ctry.lims$total>0] = 0
#  ctry.lims$total = ctry.lims$total * ((dem.lnd+spa.lnd)/spa.lnd)
#  ctry.lims$total = -ctry.lims$total
#  
  # Correcting country limits units from ha to sq.km
#  ctry.lims$total =  ctry.lims$total / 100 
#}

# Constrained-scenario (scc) suffix
suffix.list = list(unconstrained = "", country.limits = "")
                   #country.limits = "_ctrylim")

# Benchmarking scenarios (scb) #################################################
slist.names = list("scen_cb-oc", "scen_cb", "scen_bd-oc", "scen_bd",
                   "scen_cb-bd-oc", "scen_cb-bd", "scen_oc", "scen_rnd",
                   "scen_world")
slist.names = paste0(slist.names, ublim.suffix)
slist.sense = list("<=", "<=", "<=", "<=", "<=", "<=", "<=", NA, "<=")
slist.g_scalar_area = list(1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4, 1.e-4, 1, 1.e-4)
#slist.nsteps = list(1, 1, 10, 10, 10, 10, 1, 10, 10)
slist.nsteps = list(1, 1, 10, 10, 10, 10, 1, 30, 20)
slist.form = list("cb-oc", "cb", "bd-oc", "bd", "cb-bd-oc", "cb-bd", "oc",
                  "rnd", "wrld")
slist.weights = list(list(w.cb=1,w.bd=1), list(w.cb=1,w.bd=1),
                     list(w.cb=1,w.bd=1), list(w.cb=1,w.bd=1),
                     list(w.cb=c(2,rep(1,9)),
                          w.bd=c(1, 1, 3, 5, 10, 20, 35, 50, 100, 1000)),
                     list(w.cb=c(2,rep(1,9)),
                          w.bd=c(1, 1, 3, 5, 10, 20, 35, 50, 100, 1000)),
                     #
                     #                    list(w.cb=c(rep(1,5)),
                     #                         w.bd=c(1, 10, 100, 1000, 10000)),
                     #                    list(w.cb=c(rep(1,5)),
                     #                         w.bd=c(0.1, 0.5, 1.5, 10, 500)),
                     #
                     #                    list(w.cb=c(rep(1,9),0),
                     #                         w.bd=c(0, 1, 4, 10, 40, 100, 400, 1000, 4000, 1)),
                     #                    list(w.cb=c(rep(1,8),0),
                     #                         w.bd=c(0, 0.1, 0.5, 1.5, 4, 10, 50, 500, 1)),
                     #list(w.cb=1,w.bd=1), 1, list(w.cb=1,w.bd=10)
                     list(w.cb=c(2,rep(1,9)),
                          w.bd=c(1, 1, 3, 5, 10, 20, 35, 50, 100, 1000)),
                     1,
                     list(w.cb=c(2,rep(1,9)),
                          w.bd=c(1, 1, 3, 5, 10, 20, 35, 50, 100, 1000)))
#slist.w = list('NA', 'NA', 'NA', 'NA', 1:10, 1:9, 'NA', 'NA')

names(slist.names) = slist.form; names(slist.sense) = slist.form
names(slist.g_scalar_area) = slist.form; names(slist.nsteps) = slist.form
names(slist.form) = slist.form; names(slist.weights) = slist.form
#names(slist.w) = slist.form

# Control of the scenarios to be run
sct.range = names(rest.area)[target.range]
scc.range = names(constr.list)[ublim.cty.range]
scb.range = slist.form[bench.range]

# Calling run_scenarios function
#res.list = list()

# Defining scenario for debugging
#sct = sct.range[[1]]; scc = scc.range[[1]]; scb = scb.range[[5]]

start = Sys.time()
for (sct in sct.range){
  results.df = c()
  for (scc in scc.range){
    for (scb in scb.range){
      outdir = paste0("opt_results_", sct , "_v8/")
      display.dir = paste0("display_results_", sct , "_v8/")
      if (!dir.exists(outdir)) dir.create(outdir)
      if (!dir.exists(display.dir)) dir.create(display.dir)
      
      #if (!exists('wrld.form')){wrld.form=1}
      
      random = (scb == 'rnd')
      
      g_scalar_area = slist.g_scalar_area[scb][[1]]
      
      restoration.area = rest.area[[sct]]
      restoration.rhs = (restoration.area / A) * g_scalar_area
      rhs.list = list(unconstrained = restoration.rhs,
                      country.limits = c((ctry.lims$total / A) * g_scalar_area,
                                         restoration.rhs))

      slist.names = lapply(slist.names, function(x){paste0(x,suffix.list[scc])})
      constr = constr.list[scc][[1]]
      
      save.name = slist.names[scb][[1]]
      scen = save.name
      objfuncform = scb
      sense = rep(slist.sense[scb][[1]], nrow(constr))
      
      nsteps = slist.nsteps[scb][[1]]
      if ((nsteps>1) & (exists('overwrite.nsteps'))){nsteps=overwrite.nsteps}
      
      #if (exists('flat.ctrylim')) {rhs = rhs.list[scc][[1]]
      #                             rhs[length(rhs)] = rhs[length(rhs)]/nsteps
      #                            } else {rhs = rhs.list[scc][[1]] / nsteps}
      rhs = rhs.list[scc][[1]] * ifelse(random, 1, (1 / nsteps))
      
      if(exists('wgt.range')){
        weights = as.list(as.data.frame(slist.weights[scb][[1]])[wgt.range,])}
      else {weights = slist.weights[scb][[1]]}
      #if(exists("flat.ctrylim")){weights = slist.weights[5][[1]]}
      
      weightm = matrix(unlist(weights), ncol = length(weights))
      save(weightm, file=paste0(outdir, save.name, "_weightm.RData"))
      n.weights = ifelse(random, 1, length(weights$w.cb))
      
      source("run_scenario_fun.r")
      
      #res.list[scen.names[scb][[1]]] = result
      #res.list[scen.names[scb][[1]]] = run_scenario_fun(bd, cb, oc, scen.form = scb,
      #                                                  ub, constr, sense = iter.sense,
      #                                                  rhs, nsteps = iter.nsteps,
      #                                                  weights = iter.weights,
      #                                                  save.name = iter.name,
      #                                                  random = (scb == 'rnd')) 
    }
  }  
  save(results.df, file=paste0(outdir,"allscenarios_results.df_",ublim.suffix,"_",Sys.Date(), ".RData"))
  write.csv(results.df, file=paste0(display.dir,"allscenarios_results.df_",ublim.suffix,"_",Sys.Date(), ".csv"))
}



#load('./allscenarios_results.df_-_-oc_-_2018-08-15.RData')

