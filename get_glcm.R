#---------------------------------------------#
#------- Calculate Textures from GLCMs -------#
#---------------------------------------------#

#--- Goal:
#---  - Calculate texture properties from each band

#--- load libs
library(tictoc)
library(raster)
library(glcm)
library(stringr)
library(rgdal)

#--- wd
wd = paste0(getwd(),'/')

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

#--- parameters
plot_screen = F
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

#--- read dashboard of runs
dash = read.csv(paste0(wd,'jobs/glcm/dashboard_glcm.csv'), as.is = T)

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
          
          #--- read shapefile of the plots
          shp_plots = readOGR(dsn= paste0(wd,"BIomass_Samples_Shapefiles"),
                              layer=paste0("PAU_BS_",actual_dates_simple$shapefile),
                              verbose=F)
          
          #--- extract mean GLCM by plot
          message('Extracting mean GLCM by plot ...')
          GLCM_plots_mean = extract(glcm, shp_plots, fun = mean, df = T, na.rm = T, sp = T, weights = T)
          
          #--- write
          GLCM_df = GLCM_plots_mean@data
          write.csv(GLCM_df, paste0(wd,'Results/GLCM/GLCM_',ws,'_',di,'_',ba,'_',hb,'.csv'), row.names=F)
          toc()
        }
      }
    }
  }
}

message('done!')

