---

# PauWetGrassML

Developing machine learning models (RF and PLS) for predicting species richness in a wet grassland field under varying cutting systems in north Germany (Paulinenaue). All ML models and data processing was done using the [R language](https://www.r-project.org/) version 4.0.0.

---

### Rscript Workflow

1. Calculate and merge vi, dem, glcm features: get_vi.R, get_glcm.R, merge_features.R
   
2. Re-scale features for spatial analysis: get_scaled_vi.R, get_scaled_dem.R, get_scaled_glcm.R

3. Resample all spatial features to plot-scale: resample_rasters.R
   
4. Calibration, cross-validation and feature importance for RF and PLS models: ML_models_reps_tune.R
   
5. Make maps for spatially predicted species richness: make_raster_maps.R

![workflow](https://github.com/user-attachments/assets/17cc5545-7058-4c9d-98db-fd8068fbf62a)

---

### Data and Folder Structure:

BIomass_Samples_Shapefiles/:
- shapefiles of plots contours in the field

TIF/:
- geotiff files of DEM/DTM, multispectral data, and intermediate geotif files with vi,glcm,dem features.

dates.csv:
- the corresponding dates for each biomass sample/shapefile/drone flight

lib/:
- custom R functions

jobs/:
- configuration files for cluster run

Results/:
- Final and intermediate results folder
- Results/ALLDATA/: merged features 
- Results/ALLDATA/MODELPERF/: cross-validation and feature importance results
- Results/ALLDATA/MODELPERF/RASTER/: predicted spatially explicit maps of species richness

---

### Funding: 

- [DAKIS](https://adz-dakis.com/en/)

---


