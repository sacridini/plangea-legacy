rm(list=ls())

setwd('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_priorization/inputdata_v6/')

library(raster)

r.bg = raster('../rawdata/5km/others/background_frame_World_4.9km_Molweide.tif')

load('./A.RData')

load('./terrestrial_index.RData')
r.terr = r.bg
r.terr[terrestrial_index] = 1
r.terr[r.terr!=1] = NA
r.terr[r.terr==1] = 0

load('./spid_proc.RData')
load('./sphabm_proc.RData')
load('./species_index_list_terr_proc.RData')

load('./prop.hab.t0.RData')
load('./prop.hab.oa.RData')

savedir = '../rawdata/5km/BD_ranges/preprocessed/'
dir.t0 = 'current'
dir.oa = 'original'

for (i in 1:dim(sphabm_proc)[1]){
#for (i in 1:10){
  print(paste0('Habitat rasters for species ID ', spid_proc[[i]], ' under work (',
               round(100*i/dim(sphabm_proc)[1],digits=3),'% done)'))
  hab.ras.t0 = r.terr
  hab.ras.oa = r.terr
  for(j in 1:5){
    if (sphabm_proc[i,j] == 1){
      sp.ptr = terrestrial_index[species_index_list_terr_proc[[i]]]
      
      hab.vals.t0 = prop.hab.t0[species_index_list_terr_proc[[i]], j]
      hab.ras.t0[sp.ptr] = hab.ras.t0[sp.ptr] + hab.vals.t0
      
      hab.vals.oa = prop.hab.oa[species_index_list_terr_proc[[i]], j]
      hab.ras.oa[sp.ptr] = hab.ras.oa[sp.ptr] + hab.vals.oa
      
      #print(paste(length(sp.ptr), length(species_index_list_terr_proc[[i]]),
      #            length(hab.vals.t0), length(hab.vals.oa),
      #            round(mean(hab.vals.t0),digits=2), round(mean(hab.vals.oa),digits=2)))
    }
  }
  #print(spplot(stack(hab.ras.t0, hab.ras.oa), names.attr = c('Current', 'Original')))
  writeRaster(hab.ras.t0, filename=paste0(savedir,dir.t0,'/',spid_proc[[i]],'.tif'), overwrite=T)
  writeRaster(hab.ras.oa, filename=paste0(savedir,dir.oa,'/',spid_proc[[i]],'.tif'), overwrite=T)
}


# Checking consistency

test.id = spid_proc[unique(trunc(runif(200)*length(spid_proc)))]

a.id = dir('../rawdata/5km/BD_ranges/raster_amphibians/', pattern='.tif')
b.id = dir('../rawdata/5km/BD_ranges/raster_birds/', pattern='.tif')
m.id = dir('../rawdata/5km/BD_ranges/raster_mammals/', pattern='.tif')

errors = 0

for (i in test.id){
  # range raster
  if (paste0(i,'.tif') %in% a.id) {r.ras = raster(paste0('../rawdata/5km/BD_ranges/raster_amphibians/',i,'.tif'))}
  if (paste0(i,'.tif') %in% b.id) {r.ras = raster(paste0('../rawdata/5km/BD_ranges/raster_birds/',i,'.tif'))}
  if (paste0(i,'.tif') %in% m.id) {r.ras = raster(paste0('../rawdata/5km/BD_ranges/raster_mammals/',i,'.tif'))}
  # current raster
  t0.ras = raster(paste0(savedir,'current/',i,'.tif'))
  # original-area raster
  oa.ras = raster(paste0(savedir,'original/',i,'.tif'))
  # test fails if a pixel of the sum has value between 0 and 1
  if (length(which(values(((r.ras + t0.ras) > 0) & ((r.ras + t0.ras) < 1))))){errors = errors+1}
  if (length(which(values(((r.ras + oa.ras) > 0) & ((r.ras + oa.ras) < 1))))){errors = errors+1}
  print(errors)
}
