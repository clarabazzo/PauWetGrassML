RF_model = function(alldata,
                    pred_var,
                    targ_var,
                    ntree,
                    vname){
  
  if(!'train' %in% names(alldata)){
    stop('The dataframe "alldata" must contain a boolean column named "train" indicating the trainning subset (=TRUE)')
  }
  
  #--- default
  if(missing(ntree)){ntree = 100} # from rf examples
  
  #--- subset data for respective target and features
  data = alldata[,c(key_vars,targ_var,pred_var,'train')]
  
  diff_NA = nrow(na.omit(data[,c(pred_var,targ_var)])) - nrow(data[,c(pred_var,targ_var)])
  if(diff_NA != 0){
    message('NAs removed so training data has ',diff_NA,' difference in row numbers')
    data = na.omit(data[alldata$data,c(key_vars,targ_var,pred_var,'train')])
  }
  
  #--- training set
  train = data[data$train, ]
  
  #--- adjust randomForest model with all features
  rf.fit = randomForest(x = train[,pred_var],             # predictors 
                        y = train[,targ_var],             # target variable
                        data=train,                       # train dataset
                        ntree=100,                        # nr of trees in the 'forest'
                        mtry=max(1,length(pred_var)/3),   # nr of predictors(or variables) that are randomly picked to grow the threes
                        keep.forest=TRUE,                 # keep forest model object (to be used in further steps)
                        importance=TRUE)                  # assess importance of each predictor
  
  #--- Apply for the whole dataset
  data$sim_RF = predict(rf.fit, data)
  if(grepl('_obs', targ_var)){
    sim_varname = gsub('_obs','_sim',targ_var)
  }else{
    sim_varname = paste0(targ_var,'_sim')
  }
  colnames(data)[colnames(data) == 'sim_RF'] = sim_varname
  
  #--- get performances
  perf_features_train = 
    mperf(sim  = data[data$train, sim_varname],
          obs  = data[data$train, targ_var], vnam = vname, dchart = F)
  
  perf_features_test = 
    mperf(sim  = data[!data$train, sim_varname],
          obs  = data[!data$train, targ_var], vnam = vname, dchart = F)
  
  perf_features_train$subset = 'train'
  perf_features_test$subset = 'test'
  
  message(length(pred_var), ' Features, RMSE = ',round(perf_features_train$rmse,3))
  
  #--- Get variable importance from the model fit
  #--- got from here: https://hackernoon.com/random-forest-regression-in-r-code-and-interpretation
  ImpData = as.data.frame(importance(rf.fit))
  ImpData$Var.Names = row.names(ImpData)
  ImpData = ImpData[order(ImpData$`%IncMSE`, decreasing = T),]
  
  #--- bind and tag number of features and importance order
  perf_features = rbind(perf_features_train, perf_features_test)
  perf_features$total_features = length(pred_var)
  perf_features$order_features = paste(ImpData$Var.Names, collapse = '+')
  
  return(list(rf.fit = rf.fit, perf = perf_features, importance = ImpData, data = data))
  
}