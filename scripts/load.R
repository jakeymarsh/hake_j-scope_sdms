#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Environment Setup and Data Loading
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Initializes the R environment by loading required spatial and 
#'              statistical libraries, sourcing custom utility functions, 
#'              and importing processed survey and environmental data. 
#'              Defines CRS projections used for UTM Zone 8 transformations.
#' #############################################################################


# 1. Load Libraries ------------------------------------------------------------
library(dplyr)
library(mgcv) 
library(pROC)
library(sf)
library(terra)
library(geosphere)
library(ggplot2)
library(gt)
library(patchwork)
library(webshot2)
library(sgeostat)
library(viridis)
library(reshape2)
library(gt)
library(chromote)

# 2. Source Functions & Projections --------------------------------------------
source('./scripts/my_functions.R')

# WGS84 (Geographic) and UTM Zone 8 (Projected)
proj_wgs <- sf::st_crs(4326)
proj_utm <- "+proj=utm +zone=8 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"

# 3. Load Spatial Assets -------------------------------------------------------
# (Bathymetry, borders, and coastline data)
load('./data/_namerica_state_utm_df.RData')
load('./data/_bath_200_utm_df.RData')
load('./data/border.RData')
load('./data/border_chull.RData')

# 4. Load Analysis Data --------------------------------------------------------
# hake_js: Survey data; temp_js: Environmental/Temperature data
load("./data/hake_js.RData")
load("./data/temp_js.RData")


# end of script

