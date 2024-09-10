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

6. Post-processing and results' analysis: plot_VIP.R, anova_denovo.R

![workflow](https://github.com/user-attachments/assets/17cc5545-7058-4c9d-98db-fd8068fbf62a)

---

### Data and Folder Structure: ###
    .
    ├── BIomass_Samples_Shapefiles [field and plots shapefiles - zenodo]
    │   └── <Shapefiles>
    ├── Results [pre-processed data and ML results - zenodo]
    │   └── ALLDATA [pre-processed input data for RF and PLS models]
    │       └── MODELPERF [RF and PLS performance and variable importance]
    │           └── RASTER [spatially-explicit predictions]    
    │   ├── GLCM [pre-processed glcm]
    │   └── VI [pre-processed vi]           
    ├── TIF [input and scaled geotiff images - zenodo]
    |   ├── rescaled
    |   └── resampled
    ├── jobs [ancillary bash and job files for cluster runs]
    ├── lib [custom R functions]    
    ├── <Rscript files>    
    ├── date.csv [cutting dates - zenodo]
    ├── merged_obs.csv [dem and nspecies data - zenodo]
    ├── README.md
    ├── LICENSE
    └── ...

Large and binary files that are not in this repository can be found in Zenodo: https://doi.org/10.5281/zenodo.13621749
    
---

### Publications:

- Bazzo et al., (2024). Integration of UAV-Sensed Features Using Machine Learning Methods to Assess Species Richness in Wet Grassland Ecosystems. Ecological Informatics. https://doi.org/10.1016/j.ecoinf.2024.102813. Tag: [Bazzo_etal_ECOINF_2024](https://github.com/clarabazzo/PauWetGrassML/tree/Bazzo_etal_ECOINF_2024)

---

### Funding: 

- [DAKIS](https://adz-dakis.com/en/)

---

