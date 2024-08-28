#------------------------#
#--- resample rasters ---#
#------------------------#

#--- resample all raster features to same resolution/extent
library(terra)
library(sf)
library(raster)

#--- wd
wd = paste0(getwd(),'/')

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

resample_ref = raster(paste0(wd,'TIF/rescaled/VI/NDVI_19062023_mean.tif'))

if(!dir.exists(paste0(wd,'TIF/resampled'))){dir.create(paste0(wd,'TIF/resampled'))}
if(!dir.exists(paste0(wd,'TIF/resampled/GLCM'))){dir.create(paste0(wd,'TIF/resampled/GLCM'))}
if(!dir.exists(paste0(wd,'TIF/resampled/VI'))){dir.create(paste0(wd,'TIF/resampled/VI'))}
if(!dir.exists(paste0(wd,'TIF/resampled/DEM'))){dir.create(paste0(wd,'TIF/resampled/DEM'))}

glcm_properties = c("mean",
               "variance",
               "homogeneity",
               "contrast",
               "dissimilarity",
               "entropy",
               "second_moment",
               "correlation")

glcm_n_bands = 7 # ignoring correlation

l_vartype = c('DEM','VI','GLCM')
# runtime per image: 6 to 7s
# DEM = 73 images
# VI = 150 images
# GLCM = 720*7[bands]=5760 images -> 11.2h

for(v in l_vartype){
  l_files = list.files(paste0(wd,'TIF/rescaled/',v), pattern = '.tif')
  for(ras_fn in l_files){
    tictoc::tic()
        
    if(v == 'GLCM'){    
      for(b in 1:glcm_n_bands){
        ras_v = raster(paste0(wd,'TIF/rescaled/',v,'/',ras_fn), band = b)    
        ras_v = resample(ras_v, resample_ref)
        writeRaster(ras_v, 
                  paste0(wd,'TIF/resampled/',v,'/',gsub('mean.tif',paste0(glcm_properties[b],'.tif'), ras_fn)), 
                  overwrite = T)  
      }    
    }else{
      ras_v = raster(paste0(wd,'TIF/rescaled/',v,'/',ras_fn))    
      ras_v = resample(ras_v, resample_ref)
      writeRaster(ras_v, 
                  paste0(wd,'TIF/resampled/',v,'/',ras_fn), 
                  overwrite = T)    
    }
    
    tictoc::toc()
  }
}
message('done')
