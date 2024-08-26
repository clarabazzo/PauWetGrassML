#------------------------------------------------------------#
#-------------- ML models WetGrassland ----------------------#
#------------------------------------------------------------#
#--- @Goal: 
#---  Calibration, cross-validation and variable importance of 
#---  RandomForest and Partial Least Squares models to predict 
#---  Wet Grassland properties under different cut systems
#---
#--- @Contact(s): mvianna@uni-bonn.de, clarabazzo@uni-bonn.de
#------------------------------------------------------------#

#--- load libs
library(ggplot2)
library(randomForest)
library(caret)
library(tibble)
library(stringr)
library(raster)
library(tictoc)

#--- wd [root of PauWetGrassML]
wd = getwd()

glcm_properties = c("mean",
               "variance",
               "homogeneity",
               "contrast",
               "dissimilarity",
               "entropy",
               "second_moment")

glcm_bands = c("R","G","B","NIR","RE")

image_dates = c("16052022",
                "14062022",
                "02082022",
                "14092022",
                "16052023",
                "06062023",
                "19062023",
                "09082023",
                "22092023")

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

DEM_properties = c('mean','sd','q50','min','max','q25','q75','q90')

#--- raster files description
log_ras = c('### Raster properties are given in the raster names separated by "_" in the following order: ',
            '',
            '1.  VariableName: biomassobs or NSpecies',
            '2.  Model: RF or PLS',
            '3.  repetition: 1,2,3,...',
            '4.  outerkfold: 1,2,3,...',
            '5.  windowsize: 3x3, 5x5,...',
            '6.  direction: 0, 45, ...',
            '7.  classfeatures: dem, dem.vi, dem.glcm,...',
            '8.  total_features: 8, 16, 35,...',
            '9.  Treatment: all, 1, 2,...',
            '10. datafilter: alldata, lodging, non-lodging',
            '11. imagedate: 22092023,...',
            '',
            'Data information match the files MODELPERF_*.csv',
            '',
            'For example: ',
            'Nspecies_RF_1_1_15x15_135_dem.vi.glcm_59_all_alldata_22092023.tif',
            '','Results last generated in: ',as.character(Sys.time()))

#--- load custom libs
invisible(sapply(list.files(path = paste0(wd,"lib/"),full.names = T),
                 function(x) source(x)))

#--- parameters
use_seed    = T             # switch fixed random seed (T/F)
spatialize  = T             # apply models to the whole field?
seedn       = 1234          # seed number
run_id      = 'alldata'     # id run
lodg_opt    = 2             # 0 = no-lodging plots, 1 = only lodging plots, 2 = all data (including lodging)
targ_var    = 'N_species'   # variable to adjust models
n_reps      = 5             # number of repetitions
Treat       = 'all'         # treatment number
tracktime   = T

#--- hyperparameter tunning
ntrees   = c(50,100,250,500)
nodesize = c(1,5)
outer_k  = 3
inner_k  = 3

#--- run setup
dash_run_fn = 'jobs/adjust_models/dashboard_adjust.csv'
dash_run    = read.csv(paste0(wd, dash_run_fn), as.is =T)

#--- texture combinations
texture_window_list = c('3x3', '5x5', '7x7', '15x15')
texture_dir_list = c(0, 45, 90, 135)

#--- index variables
key_vars    = c('id','Year','DOY','Treatment','Cut_Number')

#--- update list based on args
args = commandArgs(trailingOnly=TRUE)
if(length(args) > 0){
  joblist = as.numeric(str_split(args[1],',')[[1]])
  texture_window_list = unique(dash_run$window_size[dash_run$job %in% joblist])
  texture_dir_list = unique(dash_run$direction[dash_run$job %in% joblist])
  train_f     = unique(dash_run$train_f[dash_run$job %in% joblist])
  run_id     = unique(dash_run$run_id[dash_run$job %in% joblist])
  lodg_opt   = unique(dash_run$lodg_opt[dash_run$job %in% joblist])
  targ_var   = unique(dash_run$targ_var[dash_run$job %in% joblist])
  Treat      = unique(dash_run$Treat[dash_run$job %in% joblist])
  
  if(length(train_f) >1){stop('Only one train_f per job allowed')}
  if(length(run_id) >1){stop('Only one run_id per job allowed')}
  if(length(lodg_opt) >1){stop('Only one lodg_opt per job allowed')}
  if(length(targ_var) >1){stop('Only one targ_var per job allowed')}
  if(length(Treat) >1){stop('Only one Treat per job allowed')}
}

#--- predictors combination
#--- dem   = canopy height features
#--- vi    = vegetation indexes features
#--- glcm  = texture features
pred_var_combinations = list('dem',
                             'vi',
                             'glcm',
                             c('dem','vi'),
                             c('dem','glcm'),
                             c('vi','glcm'),
                             c('dem','vi','glcm'))

for(tw in texture_window_list){
  for(td in texture_dir_list){
    
    #--- read all data for the combination of texture window size and direction in tw,td
    alldata = read.csv(paste0(wd,'Results/ALLDATA/ALLDATA_',tw,'_',td,'.csv'), as.is = T)
    
    #--- filter treatment
    if(Treat != 'all'){
      if(!Treat %in% c(1,2,3)){stop('Treatment code not existing: ',Treat)}
      alldata = alldata[alldata$Treatment == Treat,]
    }
    
    #--- lodging filter
    # 0 = no-lodging plots, 1 = only lodging plots, 2 = all data (including lodging)
    if(lodg_opt == 0){
      alldata = alldata[alldata$Lodging == 0,]  # no-lodging
      data_filter_nm = 'nonlodging'
    }else if(lodg_opt == 1){
      alldata = alldata[alldata$Lodging >= 1,]  # only lodging
      data_filter_nm = 'lodging'  
    }else{
      alldata = alldata # all data
      data_filter_nm = 'alldata'  
    }
    
    #--- remove NAs from targ var
    diff_NA = nrow(alldata[!is.na(alldata[,targ_var]),]) - nrow(alldata)
    if(diff_NA != 0){
      warning('Target variable data contains NA that will be omitted, total rows omitted = ', diff_NA)
      alldata = alldata[!is.na(alldata[,targ_var]),]
    }
    
    #--- check ndata points per fold
    ndata_fold = nrow(alldata) / outer_k / inner_k
    if(ndata_fold < 10){
      message('--------------------')
      message('--------------------')
      message('Less than 10 datapoints per fold! this is too low and model may not be able to adjust: \n',
              'outer_kfold = ',outer_k,
              'inner_kfold = ',inner_k,
              'number of points per inner fold = ',ndata_fold)
      message('--------------------')
      message('--------------------')
      
    }
    
    #--- run for each repetition
    feature_perf_res = list() 
    
    #--- run for each feature combinations 
    for(pvar_comb in pred_var_combinations){
      
      if(tracktime){tic(msg=pvar_comb)}
      
      #--- get the features list
      if(length(pvar_comb)>1){
        
        #--- combined features
        for(pv in pvar_comb){
          if(pv == pvar_comb[1]){
            pred_var = colnames(alldata)[grepl(paste0(pv,'_'),colnames(alldata))]
          }else{
            pred_var = c(pred_var, colnames(alldata)[grepl(paste0(pv,'_'),colnames(alldata))])
          }
        }
        
      }else{
        #--- single features
        pred_var = colnames(alldata)[grepl(paste0(pvar_comb,'_'),colnames(alldata))]
      }
      
      #--- name pvar_combination
      pvar_comb_nm = paste(pvar_comb, collapse = '.')
      
      message(' ---------------------- ')
      message(' --- Running RF/PLS --- ')
      message(' ---------------------- ')
      message(' Var: ', targ_var)
      message(' Feat: ', pvar_comb_nm)
      message(' TexWS: ', tw)
      message(' TexDi: ', td)
      
      if(spatialize){
        
        if(tracktime){tic(msg='Reading rasters')}
        
        #--- get dem, vi, glcm data
        list_ras = list()
        
        message('reading glcm rasters...')
        for(gd in image_dates){
          for(gb in glcm_bands){
            for(gi in 1:length(glcm_properties)){
              glcm_ras = raster::raster(paste0(wd,'TIF/resampled/GLCM/GLCM_',tw,'_',td,'_',gb,'_',gd,'_',glcm_properties[gi],'.tif'))
              names(glcm_ras) = paste0('glcm_',glcm_properties[gi],'_',gb,'_',tw,'_',td,'_',gd)
              list_ras[[length(list_ras)+1]] = glcm_ras
              names(list_ras)[length(list_ras)] = names(glcm_ras)
            }
          }
        }
        
        #### glcm_ras@file@nbands
        message('reading vi rasters...')
        for(gd in image_dates){
          for(VI in VI_list){
            vi_ras = raster::raster(paste0(wd,'TIF/resampled/VI/',VI,'_',gd,'_mean.tif'))
            names(vi_ras) = paste0('vi_',VI,'_',gd)
            list_ras[[length(list_ras)+1]] = vi_ras
            names(list_ras)[length(list_ras)] = names(vi_ras)
          }
        }
        
        message('reading dem rasters...')
        for(gd in image_dates){
          for(dprop in DEM_properties){
            dem_ras = raster::raster(paste0(wd,'TIF/resampled/DEM/HEIGHT_DEM_',gd,'_',dprop,'.tif'))
            names(dem_ras) = paste0('dem_',dprop,'_',gd)
            list_ras[[length(list_ras)+1]] = dem_ras
            names(list_ras)[length(list_ras)] = names(dem_ras)
          }
        }
        
        #--- check output dir
        dir.create(paste0(wd,'Results/ALLDATA/MODELPERF/RASTER/'), showWarnings = F)
        write(log_ras, paste0(wd,'Results/ALLDATA/MODELPERF/RASTER/Readme.txt'))
        
        if(tracktime){toc()}
        
      }
      
      #--- get outer fold fo ko
      if(use_seed){set.seed(seedn)}
      
      perf_ko = list()
      for(reps in 1:n_reps){
        
        #--- define outer k folds for this reps
        flds = createFolds(alldata[,targ_var], k = outer_k, list = TRUE, returnTrain = FALSE)
        
        for(ko in 1:outer_k){
          
          message('Running RF: ko:',ko, ' rep:',reps)
          
          #--- get data for outer fold
          datako = alldata[flds[[paste0('Fold',ko)]],]
          
          #--- mtry for classification or regression
          if(class(datako[, targ_var]) == 'factor'){
            mtry.rf = floor(max(1, sqrt(length(pred_var))))
          }else{
            mtry.rf = floor(max(1,length(pred_var)/3))
          }
          
          #--- mtry as only one mtry.rf
          tuneGrid_mtry = expand.grid(.mtry = c(mtry.rf))
          
          #--- training control
          ctrl = trainControl(method = "cv", # cross-validation
                              number = inner_k, # number of inner folds
                              search = 'grid', # grid search
                              savePredictions = "final") # in most cases a better summary for two class problems 
          
          #--- grid for ntree and nodesizes
          tuneGrid = expand.grid(ntrees= ntrees,
                                 nodesize = nodesize)
          
          if(tracktime){tic(msg='Tunning RF')}
          
          perf_res = list()
          for(i in 1:nrow(tuneGrid)){
            
            #message('Running RF for gridseach:', i, ' ko:',ko, ' rep:',reps)
            
            #--- RF            
            rf_model = train(reformulate(pred_var, targ_var),
                             data = datako[,c(targ_var, pred_var)],
                             method = "rf",
                             importance=TRUE,
                             trControl = ctrl,
                             tuneGrid = tuneGrid_mtry,
                             ntree = tuneGrid$ntrees[i],
                             nodesize = tuneGrid$nodesize[i])            
            
            perf_i = rf_model$results
            perf_i$ntree = tuneGrid$ntrees[i]
            perf_i$nodesize = tuneGrid$nodesize[i]
            
            mperf_i = list()
            for(mtry_ki in unique(rf_model$pred$mtry)){
              for(ki in 1:inner_k){
                mperf_ki = 
                  mperf(sim = rf_model$pred$pred[rf_model$pred$Resample == paste0('Fold',ki)],
                        obs = rf_model$pred$obs[rf_model$pred$Resample == paste0('Fold',ki)], vnam = ki, dchart = F)
                mperf_ki$mtry = mtry_ki
                if(!'a' %in% names(mperf_ki)){
                  mperf_ki$a = NA; mperf_ki$b = NA
                }
                mperf_i[[length(mperf_i)+1]] = mperf_ki
              }
            }
            mperf_i = do.call(rbind, mperf_i)
            
            #--- aggregate
            mperf_i_agg = 
            data.frame(mtry = aggregate(r2 ~ mtry, mperf_i, mean)[,1],
                       r2 = aggregate(r2 ~ mtry, mperf_i, mean)[,2],
                       rmse = aggregate(rmse ~ mtry, mperf_i, mean)[,2],
                       rrmse = aggregate(rrmse ~ mtry, mperf_i, mean)[,2],
                       bias = aggregate(bias ~ mtry, mperf_i, mean)[,2],
                       ef = aggregate(ef ~ mtry, mperf_i, mean)[,2],
                       r = aggregate(r ~ mtry, mperf_i, mean)[,2],
                       d = aggregate(d ~ mtry, mperf_i, mean)[,2],
                       cc = aggregate(cc ~ mtry, mperf_i, mean)[,2],
                       n = aggregate(n ~ mtry, mperf_i, mean)[,2])
            
            perf_i = 
            merge(perf_i,
                  mperf_i_agg,
                  by='mtry')
            
            perf_res[[i]] = perf_i
            
          }
          if(tracktime){toc()}
          
          #--- bind them all together
          perf_res_RF = do.call(rbind, perf_res)
          
          #--- order by best
          perf_res_RF = perf_res_RF[order(perf_res_RF$RMSE),]
          
          #--- get the best tunning and tag ko and reps
          perf_best_RF = perf_res_RF[1,]
          perf_best_RF$ncomp = NA
          perf_best_RF$Model = 'RF'
          
          max_ncomp = length(pred_var)
          
          
          if(spatialize){
            if(tracktime){tic(msg='Applying RF spatially for each date')}
            message('Making raster maps...')
            
            #--- get best-tuned model
            rf_model = train(reformulate(pred_var, targ_var),
                             data = datako[,c(targ_var, pred_var)],
                             method = "rf",
                             importance=TRUE,
                             trControl = ctrl,
                             tuneGrid = tuneGrid_mtry,
                             ntree = perf_best_RF$ntree[1],
                             nodesize = perf_best_RF$nodesize[1])
            
            #--- make maps for each date
            for(gd in image_dates){
              
              #--- get df and valid rows
              datasp  = ras_to_df(list_ras = list_ras, list_vars = pred_var, image_date = gd, tw = tw, td = td)
              valid_rows = as.numeric(row.names(na.omit(datasp)))
              
              #--- apply model
              datasp$Predicted = NA
              datasp$Predicted[valid_rows] = predict(rf_model, datasp[valid_rows,])
              
              #--- plausability check
              datasp$Predicted[valid_rows][datasp$Predicted[valid_rows] < 0] = 0 # cap values to zero
              datasp$Predicted[valid_rows][datasp$Predicted[valid_rows] > (max(datako[,c(targ_var)]) * 2)] = NA # remove predicted values that are 2x greater than max observations
              
              #--- return values to raster
              pred_ras = list_ras[[1]]
              values(pred_ras) = datasp$Predicted
              ras_name = paste(gsub('_','',targ_var),'RF', reps, ko, tw, td,  pvar_comb_nm, length(pred_var), Treat, data_filter_nm, gd,sep = '_')
              names(pred_ras)  = ras_name
              
              writeRaster(pred_ras, 
                          paste0(wd,'Results/ALLDATA/MODELPERF/RASTER/',ras_name,'.tif'), 
                          overwrite = T)
              
            }
            if(tracktime){toc()}
          }
          
          message('Running PLS: ko:',ko, ' rep:',reps)
          
          #--- PLS
          if(tracktime){tic(msg='Tunning PLS')}
          pls_model = 
            try({train(reformulate(pred_var, targ_var),
                            data = datako[,c(targ_var, pred_var)],
                            method = "pls",
                            scale = TRUE,
                            trControl = ctrl,
                            tuneLength = max_ncomp)}, silent = T)
          
          #--- reduce max_ncomp when model cannot be tunned to higher number of components
          while(class(pls_model)[1] == 'try-error' & max_ncomp > 1){
            options(warn = 2)
            max_ncomp = max_ncomp - 1
            pls_model = 
              try({train(reformulate(pred_var, targ_var),
                         data = datako[,c(targ_var, pred_var)],
                         method = "pls",
                         scale = TRUE,
                         trControl = ctrl,
                         tuneLength = max_ncomp)}, silent = T)
            message('max_ncomp reduced in PLS to: ', max_ncomp)
          }
          options(warn = -1)
          if(tracktime){toc()}
          
          #--- best of PLS
          perf_pls = pls_model$results[pls_model$bestTune$ncomp,]
          
          #--- calculate other indexes too for best set
          mperf_f = list()
          for(ki in 1:inner_k){
            mperf_ki = 
              mperf(sim = pls_model$pred$pred[pls_model$pred$Resample == paste0('Fold',ki)],
                    obs = pls_model$pred$obs[pls_model$pred$Resample == paste0('Fold',ki)], vnam = ki, dchart = F)
            mperf_ki$mtry = mtry_ki
            if(!'a' %in% names(mperf_ki)){
              mperf_ki$a = NA; mperf_ki$b = NA
            }
            mperf_f[[length(mperf_f)+1]] = mperf_ki
          }
          mperf_f = do.call(rbind, mperf_f)
          
          #--- add other indexes to perf_pls
          perf_pls$mtry = NA
          perf_pls$r2 = mean(mperf_f$r2)
          perf_pls$rmse = mean(mperf_f$rmse)
          perf_pls$rrmse = mean(mperf_f$rrmse)
          perf_pls$bias = mean(mperf_f$bias)
          perf_pls$ef = mean(mperf_f$ef)
          perf_pls$r = mean(mperf_f$r)
          perf_pls$d = mean(mperf_f$d)
          perf_pls$cc = mean(mperf_f$cc)
          perf_pls$n = mean(mperf_f$n)
          perf_pls$ntree = NA
          perf_pls$nodesize = NA
          perf_pls$Model = 'PLS'
          
          
          if(spatialize){
            if(tracktime){tic(msg='Applying PLS spatially for each date')}
            message('Making raster maps...')
            
            #--- make maps for each date
            for(gd in image_dates){
              
              #--- get df and valid rows
              datasp  = ras_to_df(list_ras = list_ras, list_vars = pred_var, image_date = gd, tw = tw, td = td)
              valid_rows = as.numeric(row.names(na.omit(datasp)))
              
              #--- apply model
              datasp$Predicted = NA
              datasp$Predicted[valid_rows] = predict(pls_model, datasp[valid_rows,])
              
              #--- plausability check
              datasp$Predicted[valid_rows][datasp$Predicted[valid_rows] < 0] = 0 # cap values to zero
              datasp$Predicted[valid_rows][datasp$Predicted[valid_rows] > (max(datako[,c(targ_var)]) * 2)] = NA # remove predicted values that are 2x greater than max observations
              
              #--- return values to raster
              pred_ras = list_ras[[1]]
              values(pred_ras) = datasp$Predicted
              ras_name = paste(gsub('_','',targ_var),'PLS', reps, ko, tw, td,  pvar_comb_nm, length(pred_var), Treat, data_filter_nm, gd,sep = '_')
              names(pred_ras)  = ras_name
              
              writeRaster(pred_ras, 
                          paste0(wd,'Results/ALLDATA/MODELPERF/RASTER/',ras_name,'.tif'), 
                          overwrite = T)
              
            }
            if(tracktime){toc()}
          }
          
          #--- join both and tag
          perf_best = rbind(perf_pls, perf_best_RF)
          perf_best$repetition = reps
          perf_best$outer_kfold = ko
          
          perf_ko[[length(perf_ko)+1]] = perf_best
          
        }
      }
      perf_ko = do.call(rbind, perf_ko)
      
      #--- bind and tag with tw td and store to list
      perf_ko$window_size = tw
      perf_ko$direction = td
      perf_ko$class_features = pvar_comb_nm
      perf_ko$total_features = length(pred_var)
      
      #--- store
      feature_perf_res[[length(feature_perf_res)+1]] = perf_ko
      if(tracktime){toc()}
    }
    
    feature_perf_res = do.call(rbind, feature_perf_res)
    #feature_rank_res = do.call(rbind, feature_rank_res)
    
    #--- tags
    feature_perf_res$Treatment = Treat
    feature_perf_res$targ_var = targ_var
    feature_perf_res$regression_model = feature_perf_res$Model
    feature_perf_res$data_filter = data_filter_nm
    
    #feature_rank_res$Treatment = Treat
    
    #--- write perf results and settings
    write.csv(feature_perf_res, 
              paste0(wd,'Results/ALLDATA/MODELPERF/MODELPERF_',tw,'_',td,'_',run_id,'_Treat',Treat,'_reps.csv'), row.names = F)
    
    #write.csv(feature_rank_res, 
    #          paste0(wd,'Results/ALLDATA/MODELPERF/MODELRANK_',tw,'_',td,'_',run_id,'_Treat',Treat,'_reps.csv'), row.names = F)
    
    write.csv(data.frame(setting = c('use_seed','seedn','run_id','lodg_opt'),
                         value   = c(use_seed,seedn,run_id,lodg_opt)), 
              paste0(wd,'Results/ALLDATA/MODELPERF/SETTINGS_',tw,'_',td,'_',run_id,'_Treat',Treat,'_reps.csv'), row.names = F)
    
  }
}

message('done')
