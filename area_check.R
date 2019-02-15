#rm(list=ls())
library(raster)

setwd('~/Documents/IIS_PROJECTS/global_rest_prior/global_rest_priorization/')
dir = 'inputdata_v5/'

load(paste0(dir, "A.RData"))
load(paste0(dir, "terrestrial_index.RData"))

r.bg = raster('./rawdata/5km/others/background_frame_World_4.9km_Molweide.tif')

r.plc1 <- raster("./rawdata/5km/current_LU/crop_class11_2015_4.9km_Moll.tif") + r.bg
r.plc2 <- raster("./rawdata/5km/current_LU/CulGrass_class11_2015_4.9km_Moll.tif") + r.bg
r.plc3 <- raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_forest_media_4.9km_Molweide.tif") + r.bg
r.plc4 <- raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_wetlands_media_4.9km_Molweide.tif") + r.bg
r.plc5 <- raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_desert_media_4.9km_Molweide.tif") + r.bg
r.plc6 <- raster("./rawdata/5km/current_LU/NatGrass_2015_4.9km_Moll.tif") + r.bg
r.plc7 <- raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_shrubland_media_4.9km_Molweide.tif") + r.bg
r.plc8 <- raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_ice_media_4.9km_Molweide.tif") +
  raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_water_media_4.9km_Molweide.tif") +
  raster("./rawdata/5km/current_LU/ESA_landuse_300m_2015_urban_media_4.9km_Molweide.tif") + r.bg


r.teste = r.bg; r.teste[terrestrial_index] = 1

a.df = data.frame(crop=sum(r.plc1[terrestrial_index]) * A * 100)

a.df$pasture = sum(r.plc2[terrestrial_index]) * A * 100
a.df$forest = sum(r.plc3[terrestrial_index]) * A * 100
a.df$wetland = sum(r.plc4[terrestrial_index]) * A * 100
a.df$desert = sum(r.plc5[terrestrial_index]) * A * 100
a.df$grassland = sum(r.plc6[terrestrial_index]) * A * 100
a.df$shrubland = sum(r.plc7[terrestrial_index]) * A * 100
a.df$nonrest = sum(r.plc8[terrestrial_index]) * A * 100

area.df = data.frame('vector_ha' = t(a.df))
area.df$raster_ha = 0

area.df$raster_ha[1] = sum(values(r.plc1), na.rm=T) * A * 100
area.df$raster_ha[2] = sum(values(r.plc2), na.rm=T) * A * 100
area.df$raster_ha[3] = sum(values(r.plc3), na.rm=T) * A * 100
area.df$raster_ha[4] = sum(values(r.plc4), na.rm=T) * A * 100
area.df$raster_ha[5] = sum(values(r.plc5), na.rm=T) * A * 100
area.df$raster_ha[6] = sum(r.plc6[r.plc6>0], na.rm=T) * A * 100
area.df$raster_ha[7] = sum(values(r.plc7), na.rm=T) * A * 100
area.df$raster_ha[8] = sum(values(r.plc8), na.rm=T) * A * 100

oa.for = raster("./rawdata/5km/original_LC/OA-forest.tif") + r.bg
oa.wet = raster("./rawdata/5km/original_LC/OA-wetland.tif") + r.bg
oa.des = raster("./rawdata/5km/original_LC/OA-desert.tif") + r.bg
oa.grs = raster("./rawdata/5km/original_LC/OA-grassland.tif") + r.bg
oa.shr = raster("./rawdata/5km/original_LC/OA-shrubland.tif") + r.bg
oa.iwr = raster("./rawdata/5km/original_LC/OA-ice-water-rock.tif") + r.bg

area.df$original_ha = 0

area.df$original_ha[3] = sum(values(oa.for), na.rm=T) * A * 100
area.df$original_ha[4] = sum(values(oa.wet), na.rm=T) * A * 100
area.df$original_ha[5] = sum(values(oa.des), na.rm=T) * A * 100
area.df$original_ha[6] = sum(values(oa.grs), na.rm=T) * A * 100
area.df$original_ha[7] = sum(values(oa.shr), na.rm=T) * A * 100
area.df$original_ha[8] = sum(values(oa.iwr), na.rm=T) * A * 100

area.df$vec_to_ras_diff = area.df$vector_ha - area.df$raster_ha

area.df = rbind(area.df, 'TOTAL' = colSums(area.df))

write.csv(area.df, './area_check_4.9km_Molweide3.csv')
