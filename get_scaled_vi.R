#----------------------------#
#--- calculate VegIndexes ---#
#----------------------------#

library(sf)
library(raster)
library(ggplot2)
library(stringr)
library(gridExtra)

#--- wd
wd = paste0(getwd(),'/')

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

if(!dir.exists(paste0(wd,'TIF/rescaled'))){dir.create(paste0(wd,'TIF/rescaled'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/GLCM'))){dir.create(paste0(wd,'TIF/rescaled/GLCM'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/VI'))){dir.create(paste0(wd,'TIF/rescaled/VI'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/DEM'))){dir.create(paste0(wd,'TIF/rescaled/DEM'))}

make_VI_map = function(ras, shp, TitleText){
  
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
    labs(fill="VI") + ylab(NULL) + xlab(NULL) +
    theme_bw() + theme(legend.position = 'bottom', 
                       axis.text.y=element_text(angle=-90, hjust = 0.5))
  return(ras_gg)
}

#--- parameters
plot_screen = F
plot_pdf    = F
dates_fn    = 'dates.csv'

#--- list of indexes
VI_list     = c('BNDVI',
                'CGCI',
                'CVI',
                'EVI',
                'ExG',
                'GCI',
                'GNDVI',
                'MCARI',
                'MSAVI',
                'NDRE',
                'NDVI',
                'NGI',
                'NGRDI',
                'OSAVI',
                'RDVI',
                'SR')

#--- update list based on args
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 0){
  VI_list = str_split(args[1],',')[[1]]
}

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
  actual_dates = lapply(actual_dates, FUN = as.Date)
  
  #--- format to %dd%mm%yyyy
  actual_dates_simple <- lapply(actual_dates, FUN = format, "%d%m%Y")
  
  
  
  #--- get area from shapefiles
  shp_plots <- st_read(paste0(wd,"BIomass_Samples_Shapefiles/PAU_BS_",actual_dates_simple$shapefile,".shp"),quiet=T)
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

#--- for each VI_list
#--- re-scale images to mean plot area ~ 2.5m resolution
for(VI_name in VI_list){
  
  #--- run analysis for each date
  perf_res = list()
  ext_VI_res = list()
    
  for(d in dates$h_before){
    
    message('running analysis for ',d)
    k = dates$h_before == d
    
    #--- feed dates
    actual_dates  <- list(hbefore   = as.Date(dates$h_before[k], "%d.%m.%Y"),
                          hafter    = as.Date(dates$h_after[k], "%d.%m.%Y"),
                          height    = as.Date(dates$height[k], "%d.%m.%Y"),
                          biomass   = as.Date(dates$biomass[k], "%d.%m.%Y"),
                          shapefile = as.Date(dates$shapefile[k], "%d.%m.%Y"))
    
    #--- format to %dd%mm%yyyy
    actual_dates_simple <- lapply(actual_dates, FUN = format, "%d%m%Y")
    
    #------------------------#
    #--- Read multi-bands ---#
    #------------------------#
    
    #--- read shapefile of the plots
    shp_plots <- st_read(paste0(wd,"BIomass_Samples_Shapefiles/PAU_BS_",actual_dates_simple$shapefile,".shp"),quiet=T)
    
    #--- bands separately
    R   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=1)
    G   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=2)
    B   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=3)
    RE  = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=4)
    NIR = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=5)
    
    #--- calculate vegetation index
    message('Calculating ',VI_name,' ...')
    VI = get_VegIndex(VI_name,R,G,B,RE,NIR) # comes from lib/
    
    #--- scale factor
    mean_res = mean(res(VI))
    scale_fact = as.integer((mean_area / mean_res))
    
    #--- get VI features
    message('Extracting ',VI_name,' mean...')
    mean_VI  = raster::aggregate(VI, fact = scale_fact, fun=mean)
    
    if(plot_pdf){    
      message('Writting pdf maps ...')
      
      list_gg = list()
      list_gg[[1]] = gg_area
      
      list_gg[[length(list_gg)+1]] = 
        make_VI_map(ras = VI, shp = shp_plots,
                        TitleText = paste0(VI_name,' at original resolution (',round(mean_res,4),'m): ',d))
      list_gg[[length(list_gg)+1]] = 
        make_VI_map(ras = mean_VI, shp = shp_plots,
                        TitleText = paste0('Mean ',VI_name,' at plot resolution (',round(mean(res(mean_VI)),4),'m): ',d))
      
      ggsave(
        filename = paste0(wd,"TIF/rescaled/VI/",VI_name,"_",actual_dates_simple$hbefore,".pdf"), 
        plot = marrangeGrob(list_gg, nrow=1, ncol=1), 
        width = 7, height = 3.5, dpi = 200
        )
    }
    
    #--- save rasters
    writeRaster(mean_VI, 
                paste0(wd,'TIF/rescaled/VI/',VI_name,'_',actual_dates_simple$hbefore,'_mean.tif'), 
                overwrite = T)
    
    #--- clean memory
    rm(list=c("R","G","B","RE","NIR","VI"))
    gc()
    
  }  
}

message('done')

