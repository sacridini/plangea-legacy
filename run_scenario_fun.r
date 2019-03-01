# VERSION 0.2.2 (DRAFT August 2018) CONSISTING OF 9 COMPONENTS:
# preprocessing.r, bd_process.r, wrap_optimisation.r, optimisation.r,
# wrap_scenarios.r, run_scenario_fun.r, run_optimisation.r, postprocessing.r,
# functions.r
# Contact: Hawthorne Beyer, h.beyer@uq.edu.au
# Modified by: Alvaro Iribarrem, a.iribarrem@iis-rio.org

# This script receives specific scenario-parameter values defined inside the
# main loop in wrap_scenario.r, and decides wether to call run_optimisation.r
# or to use random allocation. It also saves the results of each solver run
# in Rdata files, and maps in png format.

# Much of the code here was inherited from the run_scenario.r script of the
# previous version. It was named as it was supposed to be a function, but
# somehow R was not passing a few parameters from wrap_scenario.r into the
# function's inner environment, so for now it is simply implemented as another
# step in the processing flow.

#run_scenario_fun = function(bd, cb, oc, scen.form, ub, constr, sense, rhs, nsteps,
#                        weights, save.name, random=F){

# (!!!) Defining constants to be used in the post-processing. In that function,
# the original values of cb and oc are re-loaded, thus it is not needed to
# correct oc / cb by their scaling factors (g_scalar_oc, and g_scalar_cb) (!!!)

# Original carbon-map unit is ton/ha, defining constant to convert to ton/sq.km
const_cb = 1E2

# Converting from ton/sq.km to Mton/sq.km
const_cb = const_cb * 1E-6

# Original opportunity costs unit is USD/ha, defining constant to convert to USD/sq.km
const_oc = 1E2

# Converting from USD/sq.km to million USD/sq.km
const_oc = const_oc * 1E-6

# Area A already given in sq.km, converting to million sq.km
const_area = 1E-6

const_bd <- g_scalar_bd


s = 1
  if (!random){
    for (w in 1:n.weights){
      w.cb = unlist(weights['w.cb'])[[w]]
      w.bd = unlist(weights['w.bd'])[[w]]
      objfun = calc.objective.function(objfuncform, bd, cb, oc, w.cb, w.bd, wrld.form=wrld.form)
      plot.pu.map(objfun * (objfun<fivenum(objfun)[4]),
                  fname=paste0(outdir, scen, "_map.objfun_", w, "_t0.png"),
                  info=info, weight.print=(n.weights>1))
      
      source("run_optimisation.r")
      
      w.prt = ifelse(n.weights>1, paste0('_weighting-', w), '')
      
      # plot map of solution
      plot.pu.map(res.grad,
                  fname=paste0(display.dir, scen, "_total.restored.pu", w.prt, ".png"),
                  info=T, area.print=F, weight.print=T, localize=(nsteps==1),
                  save.raster = T)
      #plot.pu.map(res.total.restored.pu,
      #            fname=paste0(outdir, scen, "_map.total.restored.pu", w.prt, ".png"),
      #            info=info, area.print=T, weight.print=T)
      #plot.pu.map(res.total.restored.pu / ub,
      #            fname=paste0(outdir, scen, "_map.total.restored.pu.pct", w.prt, ".png"),
      #            info=info, weight.print=T)
      
      # Saving objects
      save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
      save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
      save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
      save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
      save(res.exrisk, file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
      save(res.exrisk.sd, file=paste0(outdir, scen, "_res.exrisk.sd_w_", w, ".RData"))
      
      # Species list results
      res.spp = base.spp.list
      res.spp$HabAreaNOW = habarea.t0
      res.spp$HabAreaRES = res.spp$HabAreaNOW + res.area.restored.spp[,nsteps] #res.total.restored.spp
      res.spp$HabAreaMAX = habarea.max
      res.spp$ExtRiskNOW = extinction.risk(res.spp$HabAreaNOW, res.spp$HabAreaMAX, z=0.25)
      res.spp$ExtRiskRES = extinction.risk(res.spp$HabAreaRES, res.spp$HabAreaMAX, z=0.25)
      
      write.csv(res.spp, file=paste0(display.dir, scen, "_spp.results", w.prt, ".csv"))
      
      
      # Removing objects, performing garbage collection
      rm(res.prop.restored.pu, res.total.restored.pu, res.area.restored.spp,
         res.total.restored.spp, res.exrisk)
      gc()
    } # w
    # Including results of target [sct], under constraint [scc], in the
    # benchmark scenario [scb], in the overall 'results.df' data.frame
    if(!exists("results.df")){results.df = c()}
    results.df = rbind(results.df,
                       postprocess.results(dir, outdir, scen=save.name,
                                           quad.sd=quad.sd,
                                           print.confidence=CL.prt,
                                           hasweights=((n.weights>1))))
  } else {
    niters = nsteps
    for (k in 1:niters){
      scen <- paste0(save.name, k)
      w <- 1
      nsteps <- 1
      recs <- random.allocation(restoration.rhs, ub)
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
      res.exrisk.sd <- matrix(0, nrow=ns, ncol=nsteps)
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
      res.exrisk.sd[,1] <- extinction.risk.sd(habarea.t0 + res.total.restored.spp, habarea.max, z=0.25, z.sd=z.sd)
      res.prop.restored.pu[,1] <- result$x
      res.total.restored.pu <- result$x
      # Saving objects
      save(res.prop.restored.pu, file=paste0(outdir, scen, "_res.prop.restored.pu_w_", w, ".RData"))
      save(res.total.restored.pu, file=paste0(outdir, scen, "_res.total.restored.pu_w_", w, ".RData"))
      save(res.area.restored.spp, file=paste0(outdir, scen, "_res.area.restored.spp_w_", w, ".RData"))
      save(res.total.restored.spp, file=paste0(outdir, scen, "_res.total.restored.spp_w_", w, ".RData"))
      save(res.exrisk, file=paste0(outdir, scen, "_res.exrisk_w_", w, ".RData"))
      save(res.exrisk.sd, file=paste0(outdir, scen, "_res.exrisk.sd_w_", w, ".RData"))
      # Removing objects, performing garbage collection
      rm(res.prop.restored.pu, res.total.restored.pu, res.area.restored.spp,
         res.total.restored.spp, res.exrisk)
      gc()
    } # k
      rnd.df = c()
      for (ri in 1:niters) {
        rnd.df = rbind(rnd.df,
                       postprocess.results(dir, outdir,
                                           paste0(save.name, ri),
                                           quad.sd = quad.sd,
                                           print.confidence=CL.prt,
                                           hasweights=FALSE))
      }
      rnd.mean = scen #rnd.df[1,]
      rnd.mean[,-1] = colMeans(rnd.df[,-1])
      results.df = rbind(results.df, rnd.mean)
  }
  
#print(results.df[length(results.df[,1]),])
#  return(result)
#}