#################################################################
## Project Workflow
# 1. process_JSCOPE_temps.R
#  - processes August J-SCOPE temperature forecasts, adds shelf distance, and prepares data for training and analysis
# 2. match_data.R
#  - matches hake observations to nearest JSCOPE temperature
#  - removes all points without a nearb
# 3. make_model_formulas.R
#  - makes model formulas to be used and tested in model_selection.R
# 4. model_selection.R
#  - runs LFYO model performance analysis
#  - finds model performance metrics for each model
#  - sets best model given the lowest RMSE/highest AUC for each age class
#  - saves best model and model formula
# 5. generate_output.R
#  - runs Monte Carlo simulation for delta model
#  - calculates spatially explicit CV and mean of delta output across simulations
#  - calculates regionally aggregated mean and SE of delta output in US and CAN
# 6. publish.R
#  - makes plots for publication using saved outputs in previous scripts


# load libraries
library(dplyr)
library(gt)
library(data.table)
library(ggplot2)
library(plyr)
library(viridis)
library(sf)
library(polyclip)
library(mgcv)
library(sp)
library(sf)
library(sgeostat)
# webshot::install_phantomjs()
library(webshot2)
library(webshot)
library(sfheaders)
library(pROC)
library(ggeffects)
library(patchwork)
library(reshape2)
library(gt)
library(geosphere)
library(marmap)
library(polyclip)
# library(rgdal)
library(gridExtra)
library(terra)

setwd("C:/Users/jakem/OneDrive/Desktop/GitHub/hake_j-scope_sdms")

# source function script
source('./scripts/my_functions.R')

# load files used for plotting
load('./data/_namerica_state_utm_df.RData')
load('./data/_bath_200_utm_df.RData')
load('./data/border.RData')
load('./data/border_chull.RData')

# load files for model training/predicting
load("./data/hake_js.RData")
load("./data/temp_js.RData")


proj_wgs <- sf::st_crs(4326)
proj_utm <- "+proj=utm +zone=8 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"



