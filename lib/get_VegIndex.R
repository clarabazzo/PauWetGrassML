get_VegIndex = function(NameVE,R,G,B,RE,NIR){
  
  #----------------------------------------------------------------#
  #--- Calculate vegetation indexes based on multispectral data ---#
  #----------------------------------------------------------------#
  #----
  #--- Definitiions:
  #---  - NameVE    : Vegetation Index Name [string]
  #---  - R         : Red band [raster]
  #---  - G         : Green band [raster]
  #---  - B         : Blue band [raster]
  #---  - RE        : Red Edge band [raster]
  #---  - NIR       : Near Infrared band [raster]
  #---
  #--- Equations source: https://doi.org/10.3390/rs15030639
  #--- Contact: 
  #---  - Murilo Vianna <mvianna@uni-bonn>
  #---  - Clara Bazzo <clarabazzo@uni-bonn.de> 
  #----------------------------------------------------------------#
  
  if(NameVE == 'BNDVI'){
    VE = (NIR - B) / (NIR + B)
  }
  
  if(NameVE == 'CGCI'){
    NDVI = (NIR - R) / (NIR + R)
    VE = ((NIR - RE) / (NIR + RE)) / NDVI
  }
  
  if(NameVE == 'CVI'){
    VE = NIR / G * R / G
  }
  
  if(NameVE == 'EVI'){
    VE = 2.5 * (NIR - R) / (NIR + 6 * R - 7.5 * B + 1)
  }
  
  if(NameVE == 'ExG'){
    VE = 2 * G - R - B
  }
  
  if(NameVE == 'GCI'){
    VE = (NIR / G) - 1
  }
  
  if(NameVE == 'GNDVI'){
    VE = (NIR - G) / (NIR + G)
  }
  
  if(NameVE == 'MCARI'){
    VE = (((RE - R) - 0.2) * (RE - G)) * (RE / R)
  }
  
  if(NameVE == 'MSAVI'){
    VE = (2 * NIR + 1 - sqrt((2 * NIR + 1)^2 - 8 * (NIR - R))) / 2
  }
  
  if(NameVE == 'NDRE'){
    VE = (NIR - RE) / (NIR + RE)
  }
  
  if(NameVE == 'NDVI'){
    VE = (NIR - R) / (NIR + R)
  }
  
  if(NameVE == 'NGI'){
    VE = G / (R + G + B)
  }
  
  if(NameVE == 'NGRDI'){
    VE = (G - R) / (G + R)
  }
  
  if(NameVE == 'OSAVI'){
    VE = (NIR - R) / (NIR + R + 0.16)
  }
  
  if(NameVE == 'RDVI'){
    VE = (NIR - R) / sqrt(NIR + R)
  }
  
  if(NameVE == 'SR'){
    VE = NIR / R
  }
  
  return(VE)
}


