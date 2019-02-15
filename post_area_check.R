library(raster)

# loading input layer
cb.ras = raster('./rawdata/5km/others/C_DELTA_biomass_SOC30cm_v11.tif')

load('./inputdata_v5/prop.crop.RData')
load('./inputdata_v5/prop.cultg.RData')

# Opportunity cost lists (from layers, using master_index)
occ = raster('./rawdata/5km/others/opportunity_costs_cropland_4.9km_Molweide.tif')[master_index]
ocg = raster('./rawdata/5km/others/opportunity_costs_grassland_4.9km_Molweide.tif')[master_index]

# Combined opportunity cost
oc = occ * (prop.crop / (prop.crop+prop.cultg)) + ocg * (prop.cultg / (prop.crop+prop.cultg))


load('./inputdata_v5/master_index.RData')
load('./inputdata_v5/A.RData')

res.data = read.csv('./display_results_CBD_v7/allscenarios_results.df__2018-08-29.csv')

cb.vals = cb.ras[master_index]
load('./opt_results_CBD_v7/scen_cb_res.total.restored.pu_w_1.RData')
cb.total = sum(cb.vals * res.total.restored.pu * A)
cb.data = res.data$cb.total[res.data$scenario == 'scen_cb']
cb.ratio = cb.total/cb.data

load('./opt_results_CBD_v7/scen_oc_res.total.restored.pu_w_1.RData')
oc.total = sum(oc * res.total.restored.pu * A)
cb.data = res.data$oc.total[res.data$scenario == 'scen_oc']
oc.ratio = oc.total/oc.data

oc.data = res.data$oc.total[res.data$scenario == 'scen_oc']