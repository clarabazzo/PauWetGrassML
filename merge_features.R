#---------------------------------#
#------- Merge all results -------#
#---------------------------------#

#--- load libs

#library(rgdal)
library(ggplot2)
library(stringr)
library(reshape2)

#--- wd
wd = getwd()

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

fn_merged_obs = 'merged_obs.csv'
run_VI = F # switch to T if data was not yet extracted in Results/VI
run_GLCM = F # switch to T if data was not yet extracted in Results/GLCM

H_list = c('height_dem_mean',
           'height_dem_sd',	
           'height_dem_min',	
           'height_dem_q25',	
           'height_dem_q50',	
           'height_dem_q75',	
           'height_dem_q90',	
           'height_dem_max')

#--- list of Vegetation indexes to merge
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

#--- list of textures to merge
texture_window_list = c('3x3', '5x5', '7x7', '15x15')
texture_dir_list = c(0, 45, 90, 135)
texture_band_list= c('R','G','B','RE','NIR')

#--- get heights first
all_data = read.csv(paste0(wd,fn_merged_obs), as.is = T)

#--- key variable
key_vars = c('id','Year','DOY','Treatment','Cut_Number')

#--- filter
mer_data = all_data[,c(key_vars, 'Molehill','Lodging','N_species','height_obs','biomass_obs', H_list)]

#--- extract VI data to plot scale
if(run_VI){source('get_vi.R')}

#--- read all VIs
for(VI_name in VI_list){
  
  ext_VI_res = read.csv(paste0(wd,'Results/VI/compare_data_VI_all_',VI_name,'.csv'), as.is = T)
  
  #--- store all data into single df
  if(VI_name == VI_list[1]){
    
    #--- initialize all_VIs
    all_VIs = ext_VI_res
    all_VIs$NewVI = all_VIs$VI
    colnames(all_VIs)[colnames(all_VIs) == 'NewVI'] = VI_name
    all_VIs$VI = NULL
    all_VIs$NameVI = NULL
    
  }else{
    
    all_VIs = merge(all_VIs,
                     setNames(ext_VI_res[,c('Plot_ID','hbefore_date','VI')],
                              c('Plot_ID','hbefore_date',VI_name)),
                     by = c('Plot_ID','hbefore_date'),
                     all.x = T)
  }
}
all_VIs$id = all_VIs$Plot_ID

#--- make hbefore as Year/DOY
all_VIs$hbefore_date = as.Date(all_VIs$hbefore_date)
all_VIs$Year = as.integer(format(all_VIs$hbefore_date, '%Y'))
all_VIs$DOY  = as.integer(all_VIs$hbefore_date - as.Date(paste0(all_VIs$Year,'-01-01')) + 1)


#--- merge with heights
message('nrow before merging with VIs = ', nrow(mer_data))

mer_data = 
merge(mer_data,
      all_VIs[,c(key_vars, VI_list)],
      by = key_vars,
      all.x = T)

message('nrow after merging with VIs = ', nrow((mer_data)))

#--- now merge texture but separate file per window size and direction 
#--- like in https://doi.org/10.3390/rs12162534
if(run_GLCM){source('get_glcm.R')}

for(tw in texture_window_list){
  for(td in texture_dir_list){
    
    for(b in texture_band_list){
      
      files_tw_td = list.files(paste0(wd,'Results/GLCM/'), pattern = paste0('GLCM_',tw,'_',td,'_',b,'_'), full.names = T)
      
      tw_td_b = list()
      for(fn in files_tw_td){
        fdt = read.csv(fn, as.is = T)
        date_fn = str_split(fn,pattern = '_', simplify = T)
        date_fn = gsub('.csv','',date_fn[,ncol(date_fn)])
        date_fn = as.Date(date_fn, '%d.%m.%Y')
        fdt$date = date_fn
        tw_td_b[[length(tw_td_b)+1]] = fdt
      }
      tw_td_b = do.call(rbind, tw_td_b)
      
      colnames(tw_td_b)[grepl('glcm', colnames(tw_td_b))] = 
        paste0(colnames(tw_td_b)[grepl('glcm', colnames(tw_td_b))], '_', b)
      
      if(b == texture_band_list[1]){
        tw_td = tw_td_b
      }else{
        tw_td =
        merge(tw_td,
              tw_td_b,
              by = c('id','date'))
      }
    }
    
    tw_td$Year = as.integer(format(tw_td$date, '%Y'))
    tw_td$DOY  = as.integer(tw_td$date - as.Date(paste0(tw_td$Year,'-01-01')) + 1)
  
    mer_data_t = 
    merge(mer_data,
          tw_td,
          by = c('id','Year','DOY'),
          all.x = T)
    
    mer_data_t$date = NULL
    for(VI in VI_list){
      colnames(mer_data_t)[colnames(mer_data_t) == VI] = paste0('vi_',VI)
    }
    
    for(h in H_list){
      colnames(mer_data_t)[colnames(mer_data_t) == h] = gsub('height_','',colnames(mer_data_t)[colnames(mer_data_t) == h])
    }
    
    #--- remove glcm_correlation as it leads to several errors that cannot be handled by the ML methods
    mer_data_t = mer_data_t[,seq(1,ncol(mer_data_t))[!grepl('glcm_correlation',colnames(mer_data_t))]]
    
    write.csv(mer_data_t, 
              paste0(wd,'Results/ALLDATA/ALLDATA_',tw,'_',td,'.csv'), 
              row.names = F)
    
  }
}



