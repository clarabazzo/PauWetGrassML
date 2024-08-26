#--------------------------------#
#--- extract texture features ---#
#--------------------------------#

library(sf)
library(tictoc)
library(raster)
library(glcm)
library(stringr)
library(rgdal)
library(ggplot2)

#--- wd
wd = getwd()

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

if(!dir.exists(paste0(wd,'TIF/rescaled'))){dir.create(paste0(wd,'TIF/rescaled'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/GLCM'))){dir.create(paste0(wd,'TIF/rescaled/GLCM'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/VI'))){dir.create(paste0(wd,'TIF/rescaled/VI'))}
if(!dir.exists(paste0(wd,'TIF/rescaled/DEM'))){dir.create(paste0(wd,'TIF/rescaled/DEM'))}

make_glcm_map = function(ras, shp, TitleText){
  
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
    labs(fill="GLCM") + ylab(NULL) + xlab(NULL) +
    theme_bw() + theme(legend.position = 'bottom', 
                       axis.text.y=element_text(angle=-90, hjust = 0.5))
  return(ras_gg)
}


#--- parameters
plot_screen = F
plot_pdf    = F
dates_fn    = 'dates.csv'
joblist     = NA
n_grey      = 32

#--- update list based on args
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 0){
  joblist = as.numeric(str_split(args[1],',')[[1]])
}

#--- Corresponding dates of each TIF/Shapefile/Samples
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

#--- read dashboard of runs
dash = read.csv(paste0(wd,'jobs/rescale/GLCM/dashboard_glcm.csv'), as.is = T)

#--- filter by joblist
if(!is.na(joblist)){dash = dash[dash$jobid %in% joblist,]}

#--- prepare combinations
l_ws = unique(dash$window_size)
l_di = unique(dash$direction)
l_ba = unique(dash$band)
l_hb = unique(dash$hbefore)

for(ws in l_ws){
  for(di in l_di){
    for(ba in l_ba){
      for(hb in l_hb){
        
        f_it = 
          (dash$window_size == ws) & 
          (dash$direction == di) & 
          (dash$band == ba) & 
          (dash$hbefore == hb)
        
        if(any(f_it)){
          
          if(nrow(dash[f_it,]) > 1){warning('Repeated jobids for: ',ws,' ',di,' ',ba,' ',hb,'\nUsing first occurence')}
          
          tic()
          message('running analysis for ',hb,': ',ws,' ',di,' ',ba)
          k = dates$h_before == hb
          
          #--- feed dates
          actual_dates  <- list(hbefore   = as.Date(dates$h_before[k], "%d.%m.%Y"),
                                hafter    = as.Date(dates$h_after[k], "%d.%m.%Y"),
                                height    = as.Date(dates$height[k], "%d.%m.%Y"),
                                biomass   = as.Date(dates$biomass[k], "%d.%m.%Y"),
                                shapefile = as.Date(dates$shapefile[k], "%d.%m.%Y"))
          
          #--- format to %dd%mm%yyyy
          actual_dates_simple = actual_dates
          for(i in 1:length(actual_dates)){
            actual_dates_simple[[i]] = format(actual_dates[[i]], "%d%m%Y")
          }
          
          #-----------------#
          #--- Read band ---#
          #-----------------#
          ba_n = NA
          if(ba=='R'){ba_n = 1}
          if(ba=='G'){ba_n = 2}
          if(ba=='B'){ba_n = 3}
          if(ba=='RE'){ba_n = 4}
          if(ba=='NIR'){ba_n = 5}
          
          BAND = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=ba_n)
          
          #--- window size
          ws_c = NA
          if(ws=='3x3'){ws_c = c(3,3)}
          if(ws=='5x5'){ws_c = c(5,5)}
          if(ws=='7x7'){ws_c = c(7,7)}
          if(ws=='15x15'){ws_c = c(15,15)}
          
          #--- direction
          di_c = NA
          if(di==0){di_c = c(0,1)}
          if(di==45){di_c = c(-1,1)}
          if(di==90){di_c = c(-1,0)}
          if(di==135){di_c = c(-1,-1)}
          
          #--- get GLCM
          glcm = glcm(x=BAND, n_grey=n_grey, window=ws_c, shift=di_c)
          
          #--- scale factor
          mean_res = mean(res(glcm))
          scale_fact = as.integer((mean_area / mean_res))
          
          #--- extract mean GLCM by plot
          message('Extracting mean GLCM by plot ...')
          glcm_mean = raster::aggregate(glcm, fact = scale_fact, fun=mean)
          
          #--- save rasters
          writeRaster(glcm_mean, 
                      paste0(wd,'TIF/rescaled/GLCM/GLCM_',ws,'_',di,'_',ba,'_',actual_dates_simple$hbefore,'_mean.tif'),
                      overwrite = T)
          toc()
        }
      }
    }
  }
}

message('done!')

