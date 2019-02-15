test.prior = function(area.vals, prior.vals, gain.vals, target, opt.vals, BD=F){
  
  #index.order = master_index[order(prior.vals, decreasing=T)]
  prior.order = order(prior.vals, decreasing=T)
  area.order = area.vals[prior.order]
  area.cumsum = cumsum(area.order)
  index.sel = prior.order[area.cumsum < target]  
  
  rest.vals = rep(0, length(area.vals))
  rest.vals[index.sel] = area.vals[index.sel]
  
  if (BD) {benefit = -compute.BD(rest.vals); opt.vals=-opt.vals}
  else {benefit = rest.vals * gain.vals}
  
  print(ifelse(sum(benefit) < sum(opt.vals), "PASS", "FAIL"))
  
  return(list(selection=index.sel, impact=benefit))
}

compute.BD = function(rest.vals){
  delta.vegtype.pu <- prop.restore * rest.vals
  delta.hab.spp <- rep(0, ns)
  for (i in 1:ns){
    for(j in 1:5){
      if (sphabm_proc[i,j] == 1){
        delta.hab.spp[i] <- delta.hab.spp[i] +
          sum(delta.vegtype.pu[species_index_list_proc[[i]], j], na.rm=T)
      }
    }
  }
  delta.hab.spp <- delta.hab.spp * A
  exrisk.slope <- calc.extinction.slope(habarea.t0 + delta.hab.spp,
                                        habarea.max, z=0.25)
  map.exrisk <- extinction.risk(habarea.t0 + delta.hab.spp, habarea.max, z=0.25)
  return(map.exrisk)
}
################################################################################



target = 1500000

read.dir = paste0('./opt_results_t', target/1.e4 ,'_v5/')
load(paste0(read.dir,'scen_bd_res.exrisk_w_1.RData'))
load('./inputdata_v3/exrisk.t0.RData')

load("./inputdata_v3/A.RData")

area.vals = A * (prop.crop + prop.cultg)

bd.test = test.prior(area.vals = area.vals, prior.vals = exrisk.t0,
                     gain.vals=NA, target = target, opt.vals = res.exrisk, BD=T)



load('./opt_results_t350_v5/scen_cb_res.total.restored.pu_w_1.RData')
res.cb.unc = res.total.restored.pu
hist(res.cb.unc[res.cb.unc>0])

cb.test = test.prior(area.vals = area.vals, prior.vals = cb, gain.vals=cb,
                     target = target, opt.vals = res.cb.unc * A * cb,
                     BD=F)

rest.ras = raster('../exrisk_t0.tif')
rest.ras = rest.ras/rest.ras; rest.ras[rest.ras==1]=0
rest.ras[master_index[cb.test$selection]] = area.vals[cb.test$selection]
plot(rest.ras)

opt.ras = raster('../exrisk_t0.tif')
opt.ras = opt.ras/opt.ras; opt.ras[opt.ras==1]=0
opt.ras[master_index] = res.cb.unc * A
plot(opt.ras)

plot((rest.ras>0) - (opt.ras>0))
length(which(values((rest.ras>0) - (opt.ras>0) != 0)))
opt.ras[which(values((rest.ras>0) - (opt.ras>0) != 0))]
target - sum(area.vals[cb.test$selection])



################################################################################
# 
# val.to.ras = function(val, px.list, base.ras){
#   res = base.ras
#   res[!is.na(base.ras)] = 0
#   res[px.list] = val
#   return(res)
# }
# 
# load('./opt_results_t150_v5/scen_bd_ctrylim_res.total.restored.pu_w_1.RData')
# unc.vals = res.total.restored.pu / ub
# rest.unc = val.to.ras(unc.vals, master_index, teste)
# 
# load('./opt_results_t150_v5/scen_bd_res.total.restored.pu_w_1.RData')
# cty.vals = res.total.restored.pu / ub
# rest.cty = val.to.ras(res.total.restored.pu, master_index, teste)
# 
# #spplot(stack(rest.unc, rest.cty, rest.unc-rest.cty))
#
################################################################################
# 
# load('./opt_results_t350_v5/scen_cb_res.total.restored.pu_w_1.RData')
# res.unc = res.total.restored.pu
# 
# load('./opt_results_t350_v5/scen_cb_ctrylim_res.total.restored.pu_w_1.RData')
# res.cty = res.total.restored.pu
# 
# cb
# 
# load('./inputdata_v3/master_index.RData')
# 
# br.ptr=which(world.vals==33)
# 
# sum(res.unc[br.ptr])
