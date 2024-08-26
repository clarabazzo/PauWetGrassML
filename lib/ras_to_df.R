ras_to_df = function(list_ras, list_vars, image_date, tw, td){
  
  #--- filter by variable name
  fil_ras = list()
  for(ras_var in list_vars){
    
    if(grepl('glcm_',ras_var)){
      ras_name = paste0(ras_var,'_',tw,'_',td,'_',image_date)  
    }else{
      ras_name = paste0(ras_var,'_',image_date)  
    }
    
    if(is.null(list_ras[[ras_name]])){
      message('Problem finding raster variable ',ras_var,' for day ',image_date)
    }else{
      fil_ras[[length(fil_ras)+1]] = raster::values(list_ras[[ras_name]])
      names(fil_ras)[length(fil_ras)] = ras_var
    }
  }
  return(data.frame(do.call(cbind, fil_ras)))
}