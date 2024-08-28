#--- Read predicted maps
#library(rgdal)
library(terra)
library(sf)
library(raster)
library(ggplot2)
library(stringr)


wd = paste0(getwd(),'/')

#--- example with two rasters, you may add as many as you need in the below list
list_rasters = c('Nspecies_PLS_1_1_15x15_135_dem.vi.glcm_59_all_alldata_22092023.tif',
                 'Nspecies_RF_1_1_15x15_135_dem.vi.glcm_59_all_alldata_22092023.tif')

dir.create(paste0(wd,'Results/ALLDATA/MODELPERF/MAPS/'), showWarnings = F)

#--- Corresponding dates of each TIF/Shapefile/Samples
#---    hbefore   : DEM
#---    hafter    : DTM   
dates = read.csv(paste0(wd,dates_fn), as.is = T, sep = ';')
for(cn in names(dates)){dates[,cn] = format(as.Date(dates[,cn],"%d.%m.%Y"), "%d%m%Y")}

for(ras_name in list_rasters){
  
  message('Making map for ', ras_name)
  
  #--- raster meta info
  ras_meta = setNames(data.frame(str_split(gsub('.tif','',ras_name), '_', simplify = T)),
                      c('Variable','Model','Repetition','OutKfold','WindowSize','Direction','Features','nFeatures','Treatment','Data','Date'))
  
  #--- raster and corresponding shapefile
  ras = raster(paste0(wd,'Results/ALLDATA/MODELPERF/RASTER/',ras_name))
  shp = st_read(paste0(wd,"BIomass_Samples_Shapefiles/PAU_BS_",dates$shapefile[dates$h_before == ras_meta$Date],".shp"))
  
  #--- converts raster to df for ggplot
  ras_spdf <- as(ras, "SpatialPixelsDataFrame")
  ras_df <- as.data.frame(ras_spdf)
  colnames(ras_df) <- c("value", "x", "y")
  
  date_title = as.Date(ras_meta$Date, "%d%m%Y")
  TitleText = paste0('Predicted ',ras_meta$Variable, ' (Model:',ras_meta$Model,', Features:', ras_meta$Features,', Date:',date_title,')')
  
  ras_gg = 
    ggplot() + 
    geom_raster(data=ras_df, aes(x=x, y=y, fill=value), alpha=0.8) +
    scale_fill_gradientn(colours = rev(terrain.colors(7)), limits = c(0,max(ras_df$value, na.rm=T))) +
    geom_sf(data=shp, inherit.aes = FALSE, fill=NA)+
    ggtitle(TitleText)+
    labs(fill=ras_meta$Variable) + ylab(NULL) + xlab(NULL) +
    theme_bw() + theme(legend.position = 'bottom', 
                       axis.text.y=element_text(angle=-90, hjust = 0.5))
  
  ggsave(paste0(wd,'Results/ALLDATA/MODELPERF/MAPS/',gsub('tif','png',ras_name)),
         ras_gg,
         dpi = 500,
         height = 6,
         width = 12)
  
}
