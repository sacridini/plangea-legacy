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



target = 3500000

read.dir = paste0('./opt_results_NYDC_v6/')
load(paste0(read.dir,'scen_bd_res.exrisk_w_1.RData'))
#load('./inputdata_v3/exrisk.t0.RData')

exrisk_t0 = raster('../exrisk_t0.tif')[master_index]

load("./inputdata_v3/A.RData")

area.vals = A * (prop.crop + prop.cultg)

load(paste0(read.dir, 'scen_bd_res.total.restored.pu_w_1.RData'))
res.bd.unc = res.total.restored.pu

bd.test.a = test.prior(area.vals = area.vals, prior.vals = exrisk_t0 * area.vals,
                     gain.vals=NA, target = target, opt.vals = res.exrisk, BD=T)

bd.test.s = test.prior(area.vals = area.vals, prior.vals = exrisk_t0,
                     gain.vals=NA, target = target, opt.vals = res.exrisk, BD=T)

sum(bd.test.a$impact); sum(bd.test.s$impact); sum(res.exrisk)

# Extinction risk with restorable area
bd.ras.a = raster('../exrisk_t0.tif')
bd.ras.a = bd.ras.a/bd.ras.a; bd.ras.a[bd.ras.a==1]=0
bd.ras.a[master_index[bd.test.a$selection]] = area.vals[bd.test.a$selection]
plot(bd.ras.a)

# Extinction risk without restorable area
bd.ras.s = raster('../exrisk_t0.tif')
bd.ras.s = bd.ras.s/bd.ras.s; bd.ras.s[bd.ras.s==1]=0
bd.ras.s[master_index[bd.test.s$selection]] = area.vals[bd.test.s$selection]
plot(bd.ras.s)

# Result from optimization
opt.bd = raster('../exrisk_t0.tif')
opt.bd = opt.bd/opt.bd; opt.bd[opt.bd==1]=0
opt.bd[master_index] = res.bd.unc * A
plot(opt.bd)

opt.input = raster('../exrisk_t0.tif')
opt.input = opt.input/opt.input; opt.input[opt.input==1]=0
opt.input[master_index] = bd
plot(opt.input)

plot(raster('../exrisk_t0.tif'))





load(paste0(read.dir, 'scen_cb_res.total.restored.pu_w_1.RData'))
res.cb.unc = res.total.restored.pu
hist(res.cb.unc[res.cb.unc>0])

cb.test = test.prior(area.vals = area.vals, prior.vals = cb, gain.vals=cb,
                     target = target, opt.vals = res.cb.unc * A * cb,
                     BD=F)

cb.ras = raster('../exrisk_t0.tif')
cb.ras = cb.ras/cb.ras; cb.ras[cb.ras==1]=0
cb.ras[master_index[cb.test$selection]] = area.vals[cb.test$selection]
plot(cb.ras)

opt.cb = raster('../exrisk_t0.tif')
opt.cb = opt.cb/opt.cb; opt.cb[opt.cb==1]=0
opt.cb[master_index] = res.cb.unc * A
plot(opt.cb)

#plot((cb.ras>0) - (opt.cb>0))
length(which(values((cb.ras>0) - (opt.cb>0) != 0)))
opt.cb[which(values((cb.ras>0) - (opt.cb>0) != 0))]
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
