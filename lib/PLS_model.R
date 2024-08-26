PLS_model = function(alldata,
                     pred_var,
                     targ_var,
                     repeats,
                     tuneLength,
                     vname){
  
  if(!'train' %in% names(alldata)){
    stop('The dataframe "alldata" must contain a boolean column named "train" indicating the trainning subset (=TRUE)')
  }
  
  #--- default from caret examples 
  if(missing(tuneLength)){tuneLength = min(length(pred_var),10)}  # use max 10 or n features
  if(missing(repeats)){repeats = 5}                       # five-fold CV
  
  #--- subset data for respective target and features
  data = alldata[,c(key_vars,targ_var,pred_var,'train')]
  
  diff_NA = nrow(na.omit(data[,c(pred_var,targ_var)])) - nrow(data[,c(pred_var,targ_var)])
  if(diff_NA != 0){
    message('NAs removed so training data has ',diff_NA,' difference in row numbers')
    data = na.omit(data[alldata$data,c(key_vars,targ_var,pred_var,'train')])
  }
  
  #--- training set
  train = data[data$train, ]
  
  pls.fit = 
    caret::train(reformulate(pred_var, targ_var),
                 data = train,
                 method = "pls",
                 scale = TRUE,
                 trControl = trainControl(method="repeatedcv", repeats=repeats),
                 tuneLength = tuneLength)
  
  #--- Apply for the whole dataset
  data$sim_PLS = pls.fit %>% predict(data)
  if(grepl('_obs', targ_var)){
    sim_varname = gsub('_obs','_sim',targ_var)
  }else{
    sim_varname = paste0(targ_var,'_sim')
  }
  colnames(data)[colnames(data) == 'sim_PLS'] = sim_varname
  
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
  #--- using native caret function varImp(): https://topepo.github.io/caret/variable-importance.html
  imp = caret::varImp(pls.fit)
  ImpData = imp$importance
  ImpData$Var.Names = row.names(ImpData)
  ImpData = ImpData[order(ImpData$Overall, decreasing = T),]
  
  #--- bind and tag number of features and importance order
  perf_features = rbind(perf_features_train, perf_features_test)
  perf_features$total_features = length(pred_var)
  perf_features$order_features = paste(ImpData$Var.Names, collapse = '+')
  
  return(list(pls.fit = pls.fit, perf = perf_features, importance = ImpData, data = data))
  
}
