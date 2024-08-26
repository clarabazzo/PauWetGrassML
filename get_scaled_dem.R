#-------------------------------#
#--- extract height features ---#
#-------------------------------#

library(terra)
library(sf)
library(raster)
library(ggplot2)
library(stringr)
library(gridExtra)

#--- wd
wd = getwd()

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

if(!dir.exists(paste0(wd,'TIF/rescaled'))){dir.create(paste0(wd,'TIF/rescaled'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/GLCM'))){dir.create(paste0(wd,'TIF/rescaled/GLCM'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/VI'))){dir.create(paste0(wd,'TIF/rescaled/VI'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/DEM'))){dir.create(paste0(wd,'TIF/rescaled/DEM'))}

make_height_map = function(ras, shp, TitleText){
  
  #--- converts raster to df for ggplot
  ras_spdf <- as(ras, "SpatialPixelsDataFrame")
  ras_df <- as.data.frame(ras_spdf)
  colnames(ras_df) <- c("value", "x", "y")
  
  ras_gg = 
  ggplot() + 
    geom_raster(data=ras_df, aes(x=x, y=y, fill=value), alpha=0.8) +
    scale_fill_gradientn(colours = rev(terrain.colors(7)), limits = c(0,max(ras_df$value, na.rm=T))) +
    geom_sf(data=shp, inherit.aes = FALSE, fill=NA)+
    ggtitle(TitleText)+
    labs(fill="Height (m)") + ylab(NULL) + xlab(NULL) +
    theme_bw() + theme(legend.position = 'bottom', 
                                       axis.text.y=element_text(angle=-90, hjust = 0.5))
  return(ras_gg)
}

#--- parameters
height_threshold = 0.01
plot_screen = F
plot_pdf    = T
dates_fn    = 'dates.csv'

#--- Corresponding dates of each TIF/Shapefile/Samples
#---    hbefore   : DEM
#---    hafter    : DTM   
dates = read.csv(paste0(wd,dates_fn), as.is = T, sep = ';')

#--- get average area of all plots
area_plots = list()
for(d in dates$h_before){
  
  k = dates$h_before == d
  
  #--- feed dates
  actual_dates  <- list(hbefore   = as.Date(dates$h_before[k], "%d.%m.%Y"),
                        hafter    = as.Date(dates$h_after[k], "%d.%m.%Y"),
                        height    = as.Date(dates$height[k], "%d.%m.%Y"),
                        biomass   = as.Date(dates$biomass[k], "%d.%m.%Y"),
                        shapefile = as.Date(dates$shapefile[k], "%d.%m.%Y"))
  
  #--- convert to Rdates
  actual_dates = lapply(actual_dates, \(s)as.Date(s))
  
  #--- format to %dd%mm%yyyy
  actual_dates_simple <- lapply(actual_dates, \(s)format(s,"%d%m%Y"))
  
  #--- get area from shapefiles
  shp_plots <- st_read(paste0(wd,"BIomass_Samples_Shapefiles/PAU_BS_",actual_dates_simple$shapefile,".shp"))
  area_plots[[length(area_plots)+1]] = data.frame(h_before = d,
                                                  id = shp_plots$id,
                                                  area_m2 = as.numeric(st_area(shp_plots)))
  
}
area_plots = do.call(rbind, area_plots)

mean_area = round(mean(area_plots$area_m2),2)

gg_area = 
ggplot(area_plots, aes(x=h_before, y=area_m2)) + 
  geom_hline(aes(yintercept = mean(area_m2),colour=paste(mean_area,'m2')), linetype=2) + 
  geom_boxplot(alpha=0.5, fill='grey') + 
  labs(colour="Mean")+
  ylab("Plot area (m2)") + xlab("Dates") +
  scale_color_manual(values = c('blue')) + 
  theme_bw() + theme(legend.position = 'bottom')

#--- re-scale images to mean plot area ~ 2.5m resolution
height_thres_filter = list()
for(d in dates$h_before){
  
  message('running analysis for ',d)
  k = dates$h_before == d
  
  #--- feed dates
  actual_dates  <- list(hbefore   = as.Date(dates$h_before[k], "%d.%m.%Y"),
                        hafter    = as.Date(dates$h_after[k], "%d.%m.%Y"),
                        height    = as.Date(dates$height[k], "%d.%m.%Y"),
                        biomass   = as.Date(dates$biomass[k], "%d.%m.%Y"),
                        shapefile = as.Date(dates$shapefile[k], "%d.%m.%Y"))
  
  #--- convert to Rdates
  actual_dates = lapply(actual_dates, \(s)as.Date(s))
  
  #--- format to %dd%mm%yyyy
  actual_dates_simple <- lapply(actual_dates, \(s)format(s,"%d%m%Y"))
  
  #---------------------------------#
  #--- calculate height from DEM ---#
  #---------------------------------#
  
  #--- read shapefile of the plots
  shp_plots <- st_read(paste0(wd,"BIomass_Samples_Shapefiles/PAU_BS_",actual_dates_simple$shapefile,".shp"))
  
  #--- read DEM and DTM
  h_before <- raster(paste0(wd,"TIF/PAU_DEM_",actual_dates_simple$hbefore,".tif"))
  h_after_orig <- raster(paste0(wd,"TIF/PAU_DTM_",actual_dates_simple$hafter,".tif"))
  
  #--- resample to match h_before
  h_after <- resample(h_after_orig, h_before, filename =paste0(wd,"TIF/h_after.img"), overwrite=TRUE)
  
  #--- cut corners
  h_a <- crop(h_after, extent(h_before))
  h_b <- crop(h_before, extent(h_a))
  
  #--- layer algebra
  height <- h_b - h_a
  
  #--- remove too small heights
  tooLow_height = values(height)<height_threshold
  height_thres_filter[[length(height_thres_filter)+1]] = 
    data.frame(h_before = d,
               TotalPixels = sum(!is.na(tooLow_height)),
               PercentageTooLow = sum(tooLow_height, na.rm =T) / sum(!is.na(tooLow_height)) * 100)
  message('A total of ',round(sum(tooLow_height, na.rm =T) / sum(!is.na(tooLow_height)) * 100,2),'% of pixels were removed due to rather low heights (<',height_threshold,'m) for: ',d)
  height[tooLow_height] = NA
  
  #--- scale factor
  mean_res = mean(res(height))
  scale_fact = as.integer((mean_area / mean_res))
  
  #--- get height features
  message('Extracting height features (mean, sd, min, max...)')
  mean_height  = raster::aggregate(height, fact = scale_fact, fun=mean)
  sd_height    = raster::aggregate(height, fact = scale_fact, fun=sd)
  min_height   = raster::aggregate(height, fact = scale_fact, fun=min)
  max_height   = raster::aggregate(height, fact = scale_fact, fun=max)
  q25_height   = raster::aggregate(height, fact = scale_fact, fun = function(x, na.rm){quantile(x, 0.25, na.rm=T)})
  q50_height   = raster::aggregate(height, fact = scale_fact, fun = function(x, na.rm){quantile(x, 0.50, na.rm=T)})
  q75_height   = raster::aggregate(height, fact = scale_fact, fun = function(x, na.rm){quantile(x, 0.75, na.rm=T)})
  q90_height   = raster::aggregate(height, fact = scale_fact, fun = function(x, na.rm){quantile(x, 0.90, na.rm=T)})
  
  if(plot_pdf){
    
    message('Writting pdf maps ...')
    list_gg = list()
    list_gg[[1]] = gg_area
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = height, shp = shp_plots,
                      TitleText = paste0('Height at original resolution (',round(mean_res,4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = mean_height, shp = shp_plots,
                      TitleText = paste0('Mean height at plot resolution (',round(mean(res(mean_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = sd_height, shp = shp_plots,
                      TitleText = paste0('SD height at plot resolution (',round(mean(res(sd_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = q50_height, shp = shp_plots,
                      TitleText = paste0('Median height at plot resolution (',round(mean(res(q50_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = sd_height, shp = shp_plots,
                      TitleText = paste0('Standard deviation height at plot resolution (',round(mean(res(sd_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = min_height, shp = shp_plots,
                      TitleText = paste0('Minimum height at plot resolution (',round(mean(res(min_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = max_height, shp = shp_plots,
                      TitleText = paste0('Maximum height at plot resolution (',round(mean(res(max_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = q25_height, shp = shp_plots,
                      TitleText = paste0('Quantile25% height at plot resolution (',round(mean(res(q25_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = q75_height, shp = shp_plots,
                      TitleText = paste0('Quantile75% height at plot resolution (',round(mean(res(q75_height)),4),'m)'))
    list_gg[[length(list_gg)+1]] = 
      make_height_map(ras = q90_height, shp = shp_plots,
                      TitleText = paste0('Quantile90% height at plot resolution (',round(mean(res(q90_height)),4),'m)'))
    ggsave(
      filename = paste0(wd,"TIF/rescaled/DEM/HEIGHT_",actual_dates_simple$hbefore,".pdf"), 
      plot = marrangeGrob(list_gg, nrow=1, ncol=1), 
      width = 15, height = 9
    )
  }
  
  #--- save rasters
  writeRaster(mean_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_mean.tif'), 
              overwrite = T)
  writeRaster(sd_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_sd.tif'), 
              overwrite = T)
  writeRaster(q50_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_q50.tif'), 
              overwrite = T)
  writeRaster(min_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_min.tif'), 
              overwrite = T)
  writeRaster(max_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_max.tif'), 
              overwrite = T)
  writeRaster(q25_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_q25.tif'), 
              overwrite = T)
  writeRaster(q75_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_q75.tif'), 
              overwrite = T)
  writeRaster(q90_height, 
              paste0(wd,'TIF/rescaled/DEM/HEIGHT_DEM_',actual_dates_simple$hbefore,'_q90.tif'), 
              overwrite = T)
  
}
height_thres_filter = do.call(rbind, height_thres_filter)


#--- save filtered heights for log purposes
write.csv(height_thres_filter,
          paste0(wd,'TIF/rescaled/DEM/HEIGHT_LowFilter.csv'), 
          row.names = F)

message('done')
