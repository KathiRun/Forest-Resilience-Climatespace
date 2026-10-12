# Forest-Resilience-Climatespace

-----------------------------------------------------------------------------
This repository contains the code used for the analyses in the paper 
"Diminished resilience toward the bioclimatic limits of the world's forests"
-----------------------------------------------------------------------------

The repo is organized as follows:

Code contains all the code that was used for the analysis and to generate figures and tables

### R code
-  *00_1_Create_point_mesh.R* creates a global even point mesh
-  *00_2_Extract_climate_for_mesh.R* extracts climate values for the point mesh
-  *01_1_Climate_space_PCA.R* generates the PC1 bioclimatic axis
-  *01_2_Stratified_sampling_PC1.R* evenly samples points across 20 PC1 bins per biome
-  *01_3_Add_UUID.R* assigns UUIDs to sample_points_with_composite.csv
-  *02_1_Pivot_NDVI_to_long.R* pivots the wide NDVI dataset into a long format used for the resilience and sensitivity metric
-  *02_2_Resilience_calculation.R* calculates variance and AC1 from monthly NDVI data
-  *02_3_Build_sensitivity_input_data.R* generates the data list used for the sensitivity metric pipeline
-  *02_4_Sensitivity_metric.R* calculates the sensitivity metric from monthly NDVI and climate data
-  *03_1_Combine_datasets.R* merges all dataset needed for 03_2_Climate_space_models.R
-  *03_2_Climate_space_models.R* runs the whole analysis and generates Figures 3-6, Supplementary Figures 2-12 and Supplementary Tables 3-6
-  *03_3_SIFig.1_density_plot.R* creates the density plots in SI Figure 1
-  *Fig.2_biome_map_climate_space.R* is used to create the map in Figure 2
-  *functions_sensitivity_metric.R* holds the functions used in Sensitivity_metric.R
  
	### Gee_extraction_code
- *Build_forest_mask.js* creates the forest mask
- *sample_script_parallel_compo.py* extract the climate variables and applies the forest mask
- *sample_script_parallel_climate.py* extracts the climate time series data from Google Earth Engine
- *sample_script_parallel_ndvi.py* extracts the NDVI time series data from Google Earth Engine

  
