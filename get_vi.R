#-----------------------------#
#------- Biomass model -------#
#-----------------------------#

#--- load libs
#library(terra)
#library(sf)
library(raster)
library(rgdal)
library(ggplot2)
library(stringr)

#--- wd
wd = paste0(getwd(),'/')

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

#--- parameters
plot_screen = F
dates_fn    = 'dates.csv'
biomass_fn  = 'Paulinenaue_data_2022_2023.csv'
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

#--- for each VI_list
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
    actual_dates_simple = actual_dates
    for(i in 1:length(actual_dates)){
      actual_dates_simple[[i]] = format(actual_dates[[i]], "%d%m%Y")
    }
    
    #---------------------------------#
    #--- calculate height from DEM ---#
    #---------------------------------#
    
    #--- bands separately
    R   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=1)
    G   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=2)
    B   = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=3)
    RE  = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=4)
    NIR = raster(paste0(wd,"TIF/PAU_MULTI_",actual_dates_simple$hbefore,".tif"), band=5)
    
    #--- calculate vegetation index
    message('Calculating ',VI_name,' ...')
    VI = get_VegIndex(VI_name,R,G,B,RE,NIR)
    
    #--- read shapefile of the plots
    shp_plots = readOGR(dsn= paste0(wd,"BIomass_Samples_Shapefiles"),
                        layer=paste0("PAU_BS_",actual_dates_simple$shapefile),
                        verbose=F)
    
    #--- extract mean VI by plot
    message('Extracting mean ',VI_name,' by plot ...')
    VI_plots_mean = extract(VI, shp_plots, fun = mean, df = T, na.rm = T,sp = T, weights = T)
    
    #--- adjust df
    VI_df = VI_plots_mean@data
    colnames(VI_df) = c("Plot_ID","VI")
    VI_df$NameVI = VI_name
    
    #------------------------------------------------------------#
    #--- compare measured height (RPM) vs sensed height (DEM) ---#
    #------------------------------------------------------------#
    
    #--- read data from csv
    measured  = read.csv(paste0(wd,biomass_fn), as.is = T) 
    
    #--- convert date as Rdate
    measured$Date = as.Date(paste0(measured$Year,'-01-01')) + measured$DOY - 1
    
    #--- filter date
    measured = measured[measured$Date==actual_dates$height,]
    
    #--- remova biomass NA's
    measured = measured[!is.na(measured$EDB_gm2),]
    
    #--- join data
    ext_VI  = merge(VI_df,
                    measured,
                    by=c("Plot_ID"),
                    all.y = T)
    
    #-------------------------------------#
    #--- calculate performance indexes ---#
    #-------------------------------------#
    
    #--- list to store perf results
    perf_VI_res = list()
    
    #--- calculate statistical indexes of performance (Wallach book)
    perf_VI = 
      mperf(sim = ext_VI$VI,
            obs = ext_VI$EDB_gm2, 
            vnam = paste0(VI_name ,' ~ Biomass'),
            dchart = F)
    perf_VI$treatment = 'pooled'
    
    #--- store
    perf_VI_res[[length(perf_VI_res)+1]] = perf_VI
    
    #--- performance by treatment
    for(t in unique(ext_VI$Treatment)){
      perf_VI = 
        mperf(sim = ext_VI$VI[ext_VI$Treatment == t],
              obs = ext_VI$EDB_gm2[ext_VI$Treatment == t], 
              vnam = paste0(VI_name ,' ~ Biomass'),
              dchart = F)
      perf_VI$treatment = t
      perf_VI_res[[length(perf_VI_res)+1]] = perf_VI
    }
    perf_VI_res = do.call(rbind, perf_VI_res)
    
    #--- tag date and store
    perf_VI_res$hbefore = actual_dates_simple$hbefore
    perf_res[[length(perf_res)+1]] = perf_VI_res
    
    #--- tag and store to list
    ext_VI$hbefore = actual_dates_simple$hbefore
    ext_VI_res[[length(ext_VI_res)+1]] = ext_VI
    
  }
  ext_VI_res = do.call(rbind, ext_VI_res)
  perf_res = do.call(rbind, perf_res)
  
  #--- hbefore_date
  ext_VI_res$hbefore_date = as.Date(ext_VI_res$hbefore, "%d%m%Y")
  perf_res$hbefore_date = as.Date(perf_res$hbefore, "%d%m%Y")
  
  #--- tag VI name
  perf_res$NameVI = VI_name
  
  #--- save all results
  write.csv(perf_res,
            paste0(wd,'Results/VI/performance_VI_all_',VI_name,'.csv'), 
            row.names = F)
  
  #--- save all results
  write.csv(ext_VI_res,
            paste0(wd,'Results/VI/compare_data_VI_all_',VI_name,'.csv'), 
            row.names = F)
  
}

message('done')
