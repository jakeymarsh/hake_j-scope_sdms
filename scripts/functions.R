#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Core Functions for Data Processing, Modeling, and Visualization
#' Author: Jake Marshall
#' Date: January 2026
#' Description: This script contains modular functions for spatial data 
#'              transformation, GAM training (biomass and occurrence), 
#'              Leave-Future-Out (LFO) validation, and Monte Carlo simulations.
#' #############################################################################



# --------------------------------------------------------------------------
# Data and Pre-processing 
# --------------------------------------------------------------------------


create_vars <- function(hake_t){
  
  
  hake_t$year_fac <- as.factor(hake_t$year)
  hake_t$log_wgt_total <- ifelse(hake_t$wgt_total == 0, NA, log(hake_t$wgt_total))
  hake_t$log_age2 <- ifelse(hake_t$age2 == 0, NA, log(hake_t$age2))
  hake_t$age3_4 <- hake_t$age3 + hake_t$age4
  hake_t$log_age3_4 <- ifelse(hake_t$age3_4 == 0, NA, log(hake_t$age3_4))
  hake_t$age5_20 <- rowSums(hake_t[ , paste0('age', 5:20)])
  hake_t$log_age5_20 <- ifelse(hake_t$age5_20 == 0, NA, log(hake_t$age5_20))
  
  hake_t$present_total    <- ifelse(hake_t$wgt_total > 0, 1, 0)
  hake_t$present_age2    <- ifelse(hake_t$age2 > 0, 1, 0)
  hake_t$present_age3_4  <- ifelse(hake_t$age3_4 > 0, 1, 0)
  hake_t$present_age5_20 <- ifelse(hake_t$age5_20 > 0, 1, 0)
  
  return(hake_t)
}


# Calculate UTM coords using lat (y) and lon (x)
###########################################################################
LongLatToUTM <- function(x,y,zone){
  xy <- data.frame(ID = 1:length(x), X = x, Y = y)
  coordinates(xy) <- c("X", "Y")
  proj4string(xy) <- CRS("+proj=longlat +datum=WGS84")  ## for example
  res <- spTransform(xy, CRS(paste("+proj=utm +zone=",zone," ellps=WGS84",sep='')))
  return(as.data.frame(res))
}


add_UTM_coords <- function(data, lon.var = 'lon', lat.var = 'lat', zone = 8){
  
  dat <- data
  utm <- LongLatToUTM(dat[[lon.var]],dat[[lat.var]],zone)
  dat$easting_m <- utm$X
  dat$northing_m <- utm$Y
    
  return(dat)
}



# --------------------------------------------------------------------------
# Spatial and Bathymetric 
# --------------------------------------------------------------------------

# The following function add_shelf_distance is the work of Michael Malick
add_shelf_distance <- function(data,
                               shelf,
                               proj.data      = proj_wgs,
                               proj.shelf     = proj_wgs,
                               proj.dist      = proj_utm,
                               var.shelf.dist = "shelf_distance_m",
                               lon.var        = "lon",
                               lat.var        = "lat") {
  ## Calculate distance b/w points and a bathymetric polygon
  ##
  ## This function takes as input a data.frame of point data (i.e., `data`)
  ## and a spatial polygon object (i.e., `shelf`) and calculates the distance
  ## b/w each point and the polygon boundary. Additionally, whether the point
  ## is located inside the polygon ("off") or outside ("on") is recorded. For
  ## bathymetric polygons, "off" indicates the point is off the shelf and "on"
  ## indicates the point is on the shelf.
  ##
  ## The function returns the input `data` data.frame with one additional
  ## column:
  ##   1. var.shelf.dist = distance to shelf in units of `proj.dist`
  ##
  ##      + values = off shelf (seaward of shelf)
  ##      - values = on shelf (shorward of shelf)
  ##
  ## data = data.frame to add distance to shelf to
  ## shelf = sp object of bathymetry
  ## proj.data = projection of input `data`
  ## proj.shelf = projection of input `shelf` data, expected to be `sp` class
  ## proj.dist = projection to use for calculating distance,
  ##             distances can't be calculted in lon,lat units
  ## var.shelf.dist = name of new column in `data` specifying distance to shelf
  ## lon.var = char, column name in `data` specifying longitude
  ## lat.var = char, column name in `data` specifying latitude
  
  if(!lon.var %in% names(data))
    stop("`lon.var` not found")
  if(!lat.var %in% names(data))
    stop("`lat.var` not found")
  
  dat <- data
  
  ## Create internal row indicator
  dat$ind_shelf_internal <- 1:nrow(dat)
  
  ## Project input data
  dat.ll <- subset(dat, select = c(lon.var, lat.var, "ind_shelf_internal"))
  dat.ll <- sf::st_as_sf(dat.ll, coords = c(lon.var, lat.var))
  st_crs(dat.ll) <- proj.data
  dat.proj <- sf::st_transform(dat.ll, proj.dist)
  
  ## Project shelf data
  st_crs(shelf) <- proj.shelf
  shelf.proj <- sf::st_transform(shelf, proj.dist)
  shelf.line <- st_cast(shelf.proj, "MULTILINESTRING")
  
  ## Calc distance to bathy
  dist.all <- sf::st_distance(dat.proj, shelf.line, by_element = FALSE)
  dist.all <- matrix(dist.all, nrow = nrow(dist.all))
  dist <- apply(dist.all, 1, min)
  dat[[var.shelf.dist]] <- as.vector(dist)
  
  ## Inside bathy polygon == off shelf
  shelf.off <- dat.proj[shelf.proj, ]
  dat[[var.shelf.dist]] <- ifelse(dat$ind_shelf_internal %in%
                                    shelf.off$ind_shelf_internal,
                                  dat[[var.shelf.dist]],
                                  dat[[var.shelf.dist]] * -1)
  
  ## Clean-up
  dat$ind_shelf_internal <- NULL
  
  return(dat)
}
#########################################################################


get_bathymetry <- function(data = js_temp, lat.var = 'lat', lon.var = 'lon', res = 0.5){
  
  lon1 <- min(data[[lon.var]])
  lon2 <- max(data[[lon.var]])
  lat1 <- min(data[[lat.var]])
  lat2 <- max(data[[lat.var]])
  bath <- getNOAA.bathy(lon1, lon2, 
                        lat1, lat2, resolution = res)
}


points_within_dist <- function(p1, dist, min_years) {
  #' Identify points within 
  #' p1: Dataframe with columns: easting_m, northing_m, year
  #' dist Numeric. Threshold distance in kilometers.
  #' min_years Integer. Minimum required unique years in neighborhood.
  #' return Numeric vector of indices for points to be removed.

  
  #'   ## This function takes as input two sets of x,y locations and for each point
  ## in the first set (p1), finds the set of points in the second set (p2) that
  ## are within dist km. 
  ## Variables: 
  ## p1 has columns of "easting_m", "northing_m", and "year" in that order
  ## dist is the threshold distance between points in p1 for which they are selected
  ## min_years is the minimum number of different years of nearby data required
  ##    for a point to not be excluded
  ## The function returns 
  ## a vector ind_out with values corresponding to the index of 
  ##    points in p1 that are NOT surrounded by points from different years
  ##    These points are ideally removed outside this function
  
  ind_out <- double()
  for(i in 1:nrow(p1)){

    x.sq   <- (p1[i, 1] - p1[-c(i), 1])^2
    y.sq   <- (p1[i, 2] - p1[-c(i), 2])^2
    dist.i <- sqrt(x.sq + y.sq)

    year_i <- p1[i, 3]
    ind_i <- which(dist.i < dist)

    if(length(unique(p1[ind_i, 3])) <= min_years){
      ind_out <- c(ind_out, i)
    }
  }
  
  return(ind_out)
}



# --------------------------------------------------------------------------
# Model Training 
# --------------------------------------------------------------------------

train_js_bmass_gams <- function(data, formulas = formulas.bmass.total, subset = 'wgt_total'){
  
  ########################################
  # create a suite of gams using JSCOPE data
  # data - training data set (hake_js)
  # subset - variable that needs to be positive for training (wgt_total)
  # -- i.e. wgt_total > 0 for biomass model of log_wgt_total
  # response variable is determined in formulas
  ###########################################
  
  print('training density gam!')
  
  # Train on points with positive biomass
  train_data <- data[which(data[[subset]] > 0), ]
  
  out <- list()
  for(i in 1:length(formulas)){
    out[[i]] <- mgcv::gam(formulas[[i]],
                          data = train_data,
                          method = "REML")
  }

  return(out)
}


train_js_pa_gams <- function(data, formulas = formulas.pa.total){
  
  ########################################
  # create a suite of gams using JSCOPE data
  # data - training data set (hake_js)
  # response variable is determined in formulas
  ###########################################
  
  print('training occurrence gam!')
  
  out <- list()
  for(i in 1:length(formulas)){
    out[[i]] <- mgcv::gam(formulas[[i]],
                          data = data, 
                          family = binomial(),
                          method = "REML")
  }
  
  return(out)
}



# --------------------------------------------------------------------------
# Validation and Metrics 
# --------------------------------------------------------------------------


lfo_analysis_bmass_gams <- function(data = hake_js, 
                                     target = 'log_wgt_total', 
                                     subset = 'wgt_total', 
                                     fxn = train_js_bmass_gams, 
                                     forms = formulas$bmass.total,
                                     nyears = 7){
  
  ###############################################################################
  # Trains density gams (given by forms parameter) on first n years of data starting with n = 2
  # Predicts on n+1st year and finds rmse between prediction and observed data
  # Finds mean RMSE over all years for each gam

  years <- tail(unique(data$year), nyears)
  RMSE_list <- list()
  
  #train suite of gams using first i years of data, predict the i+1th year, and calculate the RMSE
  for(i in 1:(nyears)){
    
    data_train <- data[which(data$year < years[i]), ]
    data_next <- data[which(data$year == years[i]), ]
    
    gams_bmass_i <-  fxn(data = data_train, formulas = forms, subset = subset)
    
    # Calculates RMSE between predicted biomass and observed biomass
    rmse_i <- rmse_bmass_gams(gam_list = gams_bmass_i, prediction_data = data_next, target = target)
    
    RMSE_list[[i]] <- rmse_i
  }
  RMSE_array <- do.call(cbind, RMSE_list)
  
  mean_array <- array(dim = c(nrow(RMSE_array), ncol(RMSE_array)+1))
  mean_array[, 1:ncol(mean_array)-1] <- RMSE_array
  mean_array[,ncol(mean_array)] <- rowMeans(RMSE_array)
  RMSE_array <- mean_array 
  
  return(RMSE_array)
}


lfo_analysis_pa_gams <- function(data = hake_js, 
                                  target = 'present_total', 
                                  fxn = train_js_pa_gams, 
                                  forms = formulas$pa.total,
                                  nyears = 7){
  
  ###############################################################################
  # Trains occurrence gams (given by formulas parameter) on first n years of data starting with n = 2
  # Predicts on n+1st year and finds AUC between that prediction and the observed occurrence
  # Finds mean AUC over all years for each gam
  
  years <- tail(unique(data$year), nyears)
  AUC_list <- list()
  
  #train suite of gams using first i years of data, predict the i+1th year, and calculate the AUC
  for(i in 1:(nyears)){
    
    data_train <- data[which(data$year < years[i]), ]
    data_next <- data[which(data$year == years[i]), ]
    
    gams_pa_i <-  fxn(data = data_train, formulas = forms)
    
    # Calculates AUC between predicted presence and observed presence
    auc_i <- double(length(gams_pa_i))
    for(j in 1:length(gams_pa_i)){
      
      pred_j <- predict(gams_pa_i[[j]], newdata = data_next, type = 'response', exclude = "s(year_fac)")
      
      # AUC
      roc_curve <- roc(data_next[[target]], pred_j)
      auc_i[j] <- auc(roc_curve)
    }
    
    AUC_list[[i]] <- auc_i
  }
  AUC_array <- do.call(cbind, AUC_list)
  
  mean_array <- array(dim = c(nrow(AUC_array), ncol(AUC_array)+1))
  mean_array[, 1:ncol(mean_array)-1] <- AUC_array
  mean_array[,ncol(mean_array)] <- rowMeans(AUC_array)
  AUC_array <- mean_array 
  
  return(AUC_array)
}


rmse_bmass_gams <- function(gam_list, prediction_data, target = 'log_wgt_total'){
  ############################
  # Finds the RMSE of each of a suite of gams
  # returns a vector of length length(gam_list) holding in each gam's respective rmse
  
  rmse_all <- rep(0, length(gam_list))
  for(i in 1:length(gam_list)){
    gam_i <- gam_list[[i]]
    prediction_i <- predict(gam_i, newdata = prediction_data, type = 'response', exclude = "s(year_fac)") # Chance if year_fac is not used as random effect
    rmse_i <- rmse(prediction_i, prediction_data[[target]])
    rmse_all[i] <- rmse_i
  }
  
  return(rmse_all)
}

# 
# auc_pa_gams <- function(gam_list, prediction_data = data_next, target = 'present_total'){
#   
#   auc <- double(length(gam_list))
#   for(i in 1:length(gam_list)){
#     
#     # Predict presence under future year
#     pred_i <- predict(gam_list[[i]], prediction_data, type = 'response', exclude = "s(year_fac)")
#     
#     # Find AUC metric
#     roc_curve <- roc(prediction_data[[target]], pred_i)
#     auc_i <- auc(roc_curve)
#     
#     auc[j] <- auc_i
#   }
# }


rmse <- function(true_vec, predicted_vec){
  ############################
  # find RMSE of predicted vector as compared to the true vector
  if(length(true_vec) != length(predicted_vec)){stop('Error: trained and predicted vectors not same length')}
  
  return(sqrt(mean((true_vec - predicted_vec)^2, na.rm = T)))
}


# --------------------------------------------------------------------------
# Simulation and Output
# --------------------------------------------------------------------------

sim_delta_ntimes <- function(gam_bmass, gam_pa, pred_data, n = 5000) {

  # 1. Prediction design matrices (lpmatrix)
  Xp_bmass <- predict(gam_bmass, newdata = pred_data, type = "lpmatrix")
  Xp_pa    <- predict(gam_pa, newdata = pred_data, type = "lpmatrix")

  # Zero out s(year_fac) columns to replicate exclude = "s(year_fac)"
  bmass_yr_cols <- grep("s\\(year_fac\\)", colnames(Xp_bmass))
  if (length(bmass_yr_cols) > 0) Xp_bmass[, bmass_yr_cols] <- 0

  pa_yr_cols <- grep("s\\(year_fac\\)", colnames(Xp_pa))
  if (length(pa_yr_cols) > 0) Xp_pa[, pa_yr_cols] <- 0

  # 2. Bayesian covariance matrices of model parameters
  # unconditional = TRUE incorporates uncertainty from smoothing parameter estimation
  V_bmass <- vcov(gam_bmass, unconditional = TRUE)
  V_pa    <- vcov(gam_pa, unconditional = TRUE)

  # 3. Multivariate normal draws of coefficients
  sim_beta_bmass <- MASS::mvrnorm(n = n, mu = coef(gam_bmass), Sigma = V_bmass) # [n x p_bmass]
  sim_beta_pa    <- MASS::mvrnorm(n = n, mu = coef(gam_pa), Sigma = V_pa)       # [n x p_pa]

  # 4. Project to linear predictor space: [n_cells x n]
  eta_bmass <- Xp_bmass %*% t(sim_beta_bmass)
  eta_pa    <- Xp_pa %*% t(sim_beta_pa)

  # 5. Transform to response scale
  # Lognormal bias correction using estimated residual variance (gam_bmass$sig2)
  sim_gauss <- exp(eta_bmass + (gam_bmass$sig2 / 2))

  # Binomial logit link to probability
  sim_binomial <- plogis(eta_pa)

  # 6. Combined delta output: [n_cells x n]
  combined_sim_pred <- sim_gauss * sim_binomial
  return(combined_sim_pred)
}
# 
# 
# # test
# 
# 
# 
# # 1. Old independent draws (replicated for 500 sims)
# n_check <- 500
# pred_b_se <- predict(test_gam_bmass, newdata = test_pred, type = "link", exclude = "s(year_fac)", se.fit = TRUE)
# pred_p_se <- predict(test_gam_pa, newdata = test_pred, type = "link", exclude = "s(year_fac)", se.fit = TRUE)
# 
# old_b <- array(rnorm(nrow(test_pred) * n_check, pred_b_se$fit, pred_b_se$se.fit), dim = c(nrow(test_pred), n_check))
# old_p <- array(rnorm(nrow(test_pred) * n_check, pred_p_se$fit, pred_p_se$se.fit), dim = c(nrow(test_pred), n_check))
# old_sims <- exp(old_b + (test_gam_bmass$sig2 / 2)) * plogis(old_p)
# 
# # 2. New MVN draws (from sim_delta_ntimes)
# new_sims <- sim_delta_ntimes(test_gam_bmass, test_gam_pa, test_pred, n = n_check)
# 
# # 3. Sum across all cells (total domain biomass)
# total_old <- colSums(old_sims)
# total_new <- colSums(new_sims)
# 
# # Compare the standard deviations (the true uncertainty of regional biomass)
# cat(sprintf("Old Regional SE: %.2f\n", sd(total_old)))
# cat(sprintf("New MVN Regional SE: %.2f\n", sd(total_new)))
# cat(sprintf("Ratio (New / Old): %.2fx wider uncertainty bounds\n", sd(total_new) / sd(total_old)))
# 
# 
# 
# # ==============================================================================
# # 1. SETUP TEST ENVIRONMENT
# # ==============================================================================
# library(MASS)
# library(mgcv)
# library(ggplot2)
# 
# load('./gam/best.formulas.RData')
# 
# 
# # Pick one year of prediction data and fit a quick test model pair
# test_year <- 2015
# test_pred <- temp_js[temp_js$year == test_year, ]
# test_train <- hake_js[hake_js$year < test_year, ]
# 
# # Fit test models (or use pre-existing gam_bmass and gam_pa if already in memory)
# test_gam_bmass <- gam(best.formulas$bmass.total, 
#                       data = test_train[test_train$wgt_total > 0, ], 
#                       family = gaussian(), 
#                       method = "REML")
# 
# test_gam_pa <- gam(best.formulas$pa.total, 
#                    data = test_train, 
#                    family = binomial(), 
#                    method = "REML")
# 
# # ==============================================================================
# # 2. RUN A TEST BATCH (n = 100 is plenty for testing)
# # ==============================================================================
# n_test <- 100
# sim_mvn <- sim_delta_ntimes(
#   gam_bmass = test_gam_bmass,
#   gam_pa = test_gam_pa,
#   pred_data = test_pred,
#   n = n_test
# )
# 
# # ==============================================================================
# # 3. UNIT CHECKS: Dimensions, NAs, and Range
# # ==============================================================================
# cat("--- CHECK 1: Dimensions ---\n")
# stopifnot(nrow(sim_mvn) == nrow(test_pred))
# stopifnot(ncol(sim_mvn) == n_test)
# cat("PASSED: Matrix dimensions match grid cells x simulations.\n\n")
# 
# cat("--- CHECK 2: Numerical Validity ---\n")
# stopifnot(!any(is.na(sim_mvn)))
# stopifnot(!any(is.infinite(sim_mvn)))
# stopifnot(all(sim_mvn >= 0)) # Biomass density cannot be negative
# cat("PASSED: No NAs, Inf values, or negative densities.\n\n")
# 
# cat("--- CHECK 3: Mean Surface Consistency ---\n")
# # The average across simulation iterations should closely track mgcv's native predict()
# pred_mean_manual <- rowMeans(sim_mvn)
# 
# # Compare to standard predict() with exclude
# p_pa <- predict(test_gam_pa, newdata = test_pred, type = "response", exclude = "s(year_fac)")
# p_bmass <- exp(predict(test_gam_bmass, newdata = test_pred, type = "link", exclude = "s(year_fac)") + (test_gam_bmass$sig2 / 2))
# pred_expected <- p_pa * p_bmass
# 
# mean_corr <- cor(pred_mean_manual, pred_expected)
# cat(sprintf("Correlation between simulated mean and standard predict(): %.4f\n", mean_corr))
# stopifnot(mean_corr > 0.98)
# cat("PASSED: Simulated means align with fitted model expectation.\n\n")
# 
# # ==============================================================================
# # 4. VISUAL SMOKE TEST: Checkerboard vs Smooth Field
# # ==============================================================================
# # Plot the 1st simulation iteration to visually confirm smooth spatial gradients
# test_df <- test_pred
# test_df$sim_run_1 <- sim_mvn[, 1]
# 
# ggplot(test_df, aes(x = lon, y = lat, color = sim_run_1)) +
#   geom_point(size = 1) +
#   scale_color_viridis_c(trans = "log1p") +
#   theme_minimal() +
#   labs(
#     title = "Single Realization Surface (Iteration 1)",
#     subtitle = "Verify this looks like smooth spatial fields, NOT TV static / checkerboards",
#     color = "Density"
#   )
# 
# 
# 
# 
# sim_delta_ntimes <- function(gam_bmass, gam_pa, pred_data, n = 5000){
#   ###############################
#   ## Makes n simulations of delta output given the biomass and presence gams provided.
#   ## Returns array with columns containing full prediction sets for a given simulation and
#   ## rows containing all simulated values for a given location
#   ## Variables:
#   ## gam_bmass is the biomass gam to be used in the simulation
#   ## gam_pa is the presence gam to be used in the simulation
#   ## pred_data is the prediction dataset used to make the model output
#   ## n is the number of output simulations to be created in this function
# 
#   # Gaussian model simulation
#   ###########################################################################################
#   # Make predictions on pred_data set
#   # Simulate n gaussian output prediction sets using the fit and se.fit
#   # - simulating n sets of predictions appends each set in sequence to sim_gauss - tested
#   # restructure into array - rows are each point, columns are sets of predictions - tested
#   # transform back to normal space using Eric Ward's suggested bias correction
#   pred_gauss_with_se <- predict(gam_bmass, newdata = pred_data, type = "link", exclude = "s(year_fac)", se.fit = T)
#   pred_gauss_draws <- rnorm(length(pred_gauss_with_se$fit) * n, pred_gauss_with_se$fit, pred_gauss_with_se$se.fit) # randomly selects one value from prob dist of each prediction point based on the given fit (mean) and se
#   sim_gauss_array <- array(pred_gauss_draws, dim = c(nrow(pred_data), n))
#   sim_gauss <- exp(sim_gauss_array + c((pred_gauss_with_se$se.fit^2)/2)) # bias correction
# 
#   # Binomial model simulation
#   ###########################################################################################
#   # Make predictions on pred_data set
#   # Simulate n binomial output prediction sets using fit and se.fit
#   # restructure just like above section
#   # transform binomial to probability space
#   pred_binomial_with_se <- predict(gam_pa, newdata = pred_data, type = "link", exclude = "s(year_fac)", se.fit = T)
#   pred_binom_draws <- rnorm(length(pred_binomial_with_se$fit) * n, pred_binomial_with_se$fit, pred_binomial_with_se$se.fit)
#   sim_binom_array <- array(pred_binom_draws, dim = c(nrow(pred_data), n))
#   sim_binomial <- plogis(sim_binom_array)
# 
#   # Combine predictions
#   combined_sim_pred <- sim_gauss * sim_binomial
# 
#   return(combined_sim_pred)
# }



lfo_monte_carlo_sim <- function(train_data = hake_js, 
                                pred_data = temp_js, 
                                subset = 'wgt_total', 
                                start_train = 2012,
                                formula.bmass = best.formulas$bmass.total, 
                                formula.pa = best.formulas$pa.total, 
                                n_sims = 5000){
  
  ###############################################################################
  ## Find both Spatially Explicit and Region Averaged mean predicted biomass and SE using Monte Carlo simulation 
  ## Variables: 
  ## - train_data is the data.frame that is used in model training
  ##    - has all necessary variables to train models given the passed model formulas
  ## - pred_data is the prediction data set
  ##    - has all necessary variables to predict output given the passed model formulas
  ## - start_train is the first year for which this produces an LFO prediction It then produces one for evey subsequent year
  ## - formula.bmass is the formula for the density (biomass) model which we want to use in the delta formulation. Set to the best bmass model found in model_selection.R
  ## - formula.pa is the formula for the occurrence (presence) model which we want to use in the delta formulation. Set to the best pa model found in model_selection.R
  ## - n_sims is the number of simulations to be run for each year's LFO prediction 
  ## Returns a list
  ##  - out_reg - averaged biomass over each region (total, US, CA) for each prediction year with SD and SE of averages
  ##  - out_spat_expl - spatially explicit predictions for each prediction year (locations defined by read order, each with a unique lat/lon combination)
  
  # Set temporal boundaries for training and prediction in LFO analysis
  start_train_ind <- which(unique(pred_data$year) == start_train)
  years <- sort(unique(pred_data$year)[start_train_ind:length(unique(pred_data$year))])
  nyears <- length(years)
  
  # make lists to store regional and spat expl output
  out_reg_list <- list()
  out_spat_expl_list <- list()
  
  for(i in 1:nyears){
    
    data_train <- train_data[which(train_data$year < years[i]), ]
    data_pred <- pred_data[pred_data$year == years[i], ]
    
    gam_bmass_i <- gam(formula.bmass, 
                       data = data_train[which(data_train[[subset]] > 0), ],
                       method = "REML")
    
    gam_pa_i <- gam(formula.pa, 
                    data = data_train,
                    family = binomial(),
                    method = "REML")
    
    # simulate output n times and find mean output for each location
    simulations_i <- sim_delta_ntimes(gam_bmass = gam_bmass_i, 
                                      gam_pa = gam_pa_i, 
                                      pred_data = data_pred, 
                                      n = n_sims)
    mean_pred_i <- rowMeans(simulations_i)
    
    # find coefficient of variation (CV) at each location
    # over each row, k (each prediction location), find SE and calculate CV (se/mean)
    #################################################################################
    se_i <- double(length = nrow(simulations_i))
    cv_i <- double(length = nrow(simulations_i))
    for(k in 1:nrow(simulations_i)){
      se_i[k] <- sd(simulations_i[k, ])
      cv_i[k] <- se_i[k] / mean_pred_i[k]
    }
    out_spat_expl_list[[i]] <- data.frame(mean = mean_pred_i, 
                                          se = se_i,
                                          cv = cv_i, 
                                          year = years[i], 
                                          read_order = data_pred$read_order)
    
    
    # find SE of regionally aggregated predicted biomass
    #################################################################################
    # 1. separate predictions in US and Canada
    # 2. find mean of predictions over each region for each simulation
    # 3. find mean and sd of the region means across simulations
    
    # - predictions are made with data_pred -> output in total_sim_array corresponds to locations of pred_data
    in_US <- in.chull(data_pred$lon, data_pred$lat, border_chull$lon, border_chull$lat)
    
    # 1. Separate US and CAN simulated biomass
    US_array <- simulations_i[in_US,] 
    CAN_array <- simulations_i[!in_US,]
    
    # 2. Sum biomass across each region for each year
    US_sums <- colSums(US_array)
    CAN_sums <- colSums(CAN_array)
    
    # 3. Find how biomass sums vary across simulations
    mean_reg <- c(mean(US_sums), mean(CAN_sums))
    se_reg <- c(sd(US_sums), sd(CAN_sums))
    cv_reg <- se_reg/mean_reg
    region <- c('US', 'CAN')
    out_reg_list[[i]] <- data.frame(mean_reg, se_reg, cv_reg, region, year = years[i])
  }
  
  out_spat_expl <- do.call(rbind, out_spat_expl_list)
  out_reg <- do.call(rbind, out_reg_list)
  out_list <- list(out_spat_expl, out_reg)
  names(out_list) <- c("out_spat_expl", "out_reg")
  
  return(out_list)
}



# --------------------------------------------------------------------------
# Visualization 
# --------------------------------------------------------------------------


plot_heatmap <- function(data, years, age = 'total', names, 
                         type = 'RMSE', lims = c(0,4), data_type = 'J-SCOPE', 
                         xlab = "year predicted"){
  ############################
  # plots accuracy of each tested model in model_selection.R for each year investigated
  # RMSE is plotted for density models
  # AUC is plotted for occurrence models
  
  # flip rows for plotting
  dat <- data[rev(1:nrow(data)),]
  
  dat_df <- melt(dat)
  dat_df$model <- factor(dat_df$Var1)
  dat_df$year <- factor(dat_df$Var2)
  
  if(type == 'RMSE'){
    g <- ggplot(data = dat_df, aes(x = year, y = model, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 3)), size = 3.5) +
      scale_fill_gradientn(colours = rev(heat.colors(20)),  
                           name = expression(atop('RMSE of', paste('log'[e], '(B)'))),
                           limits = lims, 
                           oob = scales::squish) +
      scale_x_discrete(labels = years) +  
      scale_y_discrete(labels = names)  +    
      theme(text = element_text(family = "serif"),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = age,
           x = xlab,
           y = "model")
  }else{
    g <- ggplot(data = dat_df, aes(x = year, y = model, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 3)), size = 3.5) +
      scale_fill_gradientn(colours = heat.colors(20),  # You can customize the color palette
                           name = "AUC", 
                           limits = lims, 
                           oob = scales::squish) +
      scale_x_discrete(labels = years) +  # Customize X-axis labels
      scale_y_discrete(labels = names) +  # Customize Y-axis labels
      theme(text = element_text(family = "serif"),
            axis.text.x=element_text(size=10),
            axis.text.y=element_text(size=10),
            axis.title.x=element_text(size=12),
            axis.title.y=element_text(size=12),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = age,
           x = xlab,
           y = "model")
  }
  
  return(g)
}





plot_hist_years_delta <- function(data, delta_lims, years, age){
  
  ############################
  ## plot specific years' delta predictions for each age class
  ## data should be of the form spat_expl_simulation_outputs--created in generate_output.R
  
  plot_list <- list()
  
  for(i in 1:length(years)){
    plot_delta_pred <- ggplot(data[data$age == age & data$year == years[i], ]) +
      geom_point(aes(x = easting_m, y = northing_m, color = mean), shape = 15, size = 0.06) +
      geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
      geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
      labs(x = "Longitude",
           y = "Latitude",
           color = 'kg',
           title = paste0('Predicted Biomass for ', years[i])) +
      scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
      ylim(ymin, ymax) +
      scale_color_viridis(limits = delta_lims, oob = scales::squish) +
      theme_simple(base_size = 8, grid = TRUE) +
      theme(text = element_text(family = "serif"),
            plot.title = element_text(hjust = 0.5),
            panel.spacing.x = unit(-0.5, "pt"),
            panel.spacing = unit(1, "lines"),
            legend.justification = c(0, 0),
            legend.position = c(0.01, 0.005),
            legend.title.align = 0.5,
            legend.key.size = unit(0.4, "cm"),
            legend.background = element_blank(),
            legend.text = element_text(size = 6),
            legend.title = element_text(size = 6),
            panel.grid.major = element_line(size = 0.2))
    plot_list[[i]] <- plot_delta_pred
  }
  
  return(plot_list)
}




plot_ageclass_maps <- function(data, mean_lims, se_lims, change_lims, age_group){
  
  #############################
  ## Plots mean forecasted biomass maps, corresponding CV, and change in biomass from 2023 to 2024
  ## data should be of form spatExpl_sim_out--created in generate_output.R
  plot_mean <- ggplot(data[data$age == age_group, ]) +
    geom_point(aes(x = easting_m, y = northing_m, color = mean_post/1000), shape = 15, size = 0.06) +
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0('Predicted Biomass')) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_viridis(limits = mean_lims, oob = scales::squish) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(text = element_text(family = "serif"),
          plot.title = element_text(hjust = 0.5),
          panel.spacing.x = unit(-0.5, "pt"),
          panel.spacing = unit(1, "lines"),
          legend.justification = c(0, 0),
          legend.position = c(0.01, 0.005),
          legend.title.align = 0.5,
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          panel.grid.major = element_line(size = 0.2))
  
  plot_se <- ggplot(data[data$age == age_group, ]) +
    geom_point(aes(x = easting_m, y = northing_m, color = se_post/1000), shape = 15, size = 0.06) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0('Prediction SE')) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_viridis(limits = se_lims, oob = scales::squish) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(text = element_text(family = "serif"),
          plot.title = element_text(hjust = 0.5),
          panel.spacing.x = unit(-0.5, "pt"),
          panel.spacing = unit(1, "lines"),
          legend.justification = c(0, 0),
          legend.position = c(0.01, 0.005),
          legend.title.align = 0.5,
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          panel.grid.major = element_line(size = 0.2))

  
  plot_mean_change <- ggplot(data[data$age == age_group, ]) +
    geom_point(aes(x = easting_m, y = northing_m, color = mean_diff/1e3), 
               shape = 15, size = 0.06) +
    
    # Updated shaded grid logic for 25% threshold
    geom_point(data = {
      # 1. Take all points for this age group (not just significant ones)
      grid_df <- data[data$age == age_group, ]
      
      # 2. Round coordinates to create the 10km bins
      grid_df$east_grid <- round(grid_df$easting_m / 10000) * 10000
      grid_df$north_grid <- round(grid_df$northing_m / 10000) * 10000
      
      # 3. Group by the grid bins and calculate the ratio of significance
      grid_summary <- grid_df %>%
        group_by(east_grid, north_grid) %>%
        summarize(sig_ratio = mean(mean_is_significant, na.rm = TRUE), .groups = 'drop') %>%
        filter(sig_ratio >= 0.25) # The 25% threshold
      
      grid_summary
    },
    aes(x = east_grid, y = north_grid),
    shape = 3,          # Back to the vertical plus sign
    size = 0.5,         # Dropped from 0.5 to make them less crowded
    stroke = 0.2,      # Dropped from 0.3: thin enough to be "fine," thick enough to see
    color = "black",  
    alpha = 0.8) +      # Dropped from 0.5: lets the biomass colors breathe more
    
    geom_sf(data = bath_200_utm_df, fill = NA, color = "black", size = 0.2) +  
    geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey10", size = 0.2) +
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0('Change in Biomass')) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_gradient2(low = "steelblue",
                          mid = "grey90",
                          high = "red3",
                          na.value = "transparent",
                          breaks = c(min(change_lims), min(change_lims)/2, 0, 
                                     max(change_lims)/2, max(change_lims)),
                          limits = change_lims,
                          midpoint = 0,
                          oob = scales::squish) +    
    theme_simple(base_size = 8, grid = TRUE) +
    theme(text = element_text(family = "serif"),
          plot.title = element_text(hjust = 0.5),
          panel.spacing.x = unit(-0.5, "pt"),
          panel.spacing = unit(1, "lines"),
          legend.justification = c(0, 0),
          legend.position = c(0.01, 0.005),
          legend.title.align = 0.5,
          legend.key.size = unit(0.4, "cm"),
          legend.background = element_blank(),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 6),
          panel.grid.major = element_line(size = 0.2))
  
  full_plot <- gridExtra::grid.arrange(plot_mean, plot_se, plot_mean_change, ncol = 3)
  return(full_plot)
}


## The following function (theme_simple) belongs to Michael Malick
##################################################################
theme_simple <- function(base_size = 11,
                         base_family = "",
                         grid = FALSE,
                         mode = "light") {
  
  half_line <- base_size / 2
  if(mode == "dark") {
    grey0 <- "grey70"
    grey1 <- "grey30"
    grey2 <- "grey18"
    bg    <- "black"
  } else {
    grey0 <- "grey30"
    grey1 <- "grey70"
    grey2 <- "grey87"
    bg    <- NA
  }
  
  x <- theme_light(base_size = base_size, base_family = base_family) +
    theme(text = element_text(family = "serif"),
          plot.background = element_rect(fill = bg, colour = NA),
          panel.background = element_rect(fill = bg, colour = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks.length = unit(half_line / 2.2, "pt"),
          strip.background = element_rect(fill = NA, colour = NA),
          strip.text.x = element_text(colour = grey0),
          strip.text.y = element_text(colour = grey0),
          axis.text = element_text(colour = grey0),
          axis.title = element_text(colour = grey0),
          legend.title = element_text(colour = grey0, size = rel(0.9)),
          panel.border = element_rect(fill = NA, colour = grey1, size = rel(1.1)),
          legend.key.size = unit(0.9, "lines"),
          legend.text = element_text(size = rel(0.7), colour = grey0),
          legend.key = element_rect(colour = NA, fill = bg),
          legend.background = element_rect(colour = NA, fill = bg),
          plot.title = element_text(colour = grey0, size = rel(1)),
          plot.subtitle = element_text(colour = grey0, size = rel(.85)),
          plot.caption = element_text(colour = grey0, size = rel(0.8)),
          plot.tag = element_text(colour = grey0, size = rel(1.1)))
  if(grid) {
    x <- x + theme(panel.grid.major = element_line(color = grey2, size = 0.2))
  }
  return(x)
}


