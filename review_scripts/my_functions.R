
create_vars <- function(hake_t){
  
  ## Create variables in match_data.R
  hake_t$year_fac <- as.factor(hake_t$year)
  hake_t$log_wgt_total <- ifelse(hake_t$wgt_total == 0, NA, log(hake_t$wgt_total))
  hake_t$log_age2 <- ifelse(hake_t$age2 == 0, NA, log(hake_t$age2))
  hake_t$age3_4 <- hake_t$age3 + hake_t$age4
  hake_t$log_age3_4 <- ifelse(hake_t$age3_4 == 0, NA, log(hake_t$age3_4))
  hake_t$age5_20 <- rowSums(hake_t[ , paste0('age', 5:20)])
  hake_t$log_age5_20 <- ifelse(hake_t$age5_20 == 0, NA, log(hake_t$age5_20))
  
  ## Age-specific presence/absence variable
  hake_t$present_total    <- ifelse(hake_t$wgt_total > 0, 1, 0)
  hake_t$present_age2    <- ifelse(hake_t$age2 > 0, 1, 0)
  hake_t$present_age3_4  <- ifelse(hake_t$age3_4 > 0, 1, 0)
  hake_t$present_age5_20 <- ifelse(hake_t$age5_20 > 0, 1, 0)
  
  return(hake_t)
}


get_bathymetry <- function(data = js_temp, lat.var = 'lat', lon.var = 'lon', res = 0.5){
  
  lon1 <- min(data[[lon.var]])
  lon2 <- max(data[[lon.var]])
  lat1 <- min(data[[lat.var]])
  lat2 <- max(data[[lat.var]])
  bath <- getNOAA.bathy(lon1, lon2, 
                        lat1, lat2, resolution = res)
}


# bath <- get_bathymetry(data = js_temp, lat.var = 'lat_rho', lon.var = 'lon_rho', res = 0.5)
# bath.df <- bathymetry_to_df(data = bath)
# bath.df <- add_UTM_coords(data = bath.df, zone = 8)
# save(bath.df, file = './outputs/bath_df.RData')
# or load bath


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



points_within_dist <- function(p1, dist, min_years) {
  ## Find if points in p1 have another point within dist km that has another year
  ##
  ## This function takes as input two sets of x,y locations and for each point
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
    
  # find dist between p1[i] and all other points in p1 for all i
  # select those less than dist away
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

# 
#   # another method -- more efficient, less readable
#   for(i in 1:nrow(p1)){
# 
#     print(paste0(round(100*(i/nrow(p1)))))
# 
#     year_i <- p1[i, 3]
#     p2 <- p1[-c(i),]
#     years <- double()
# 
#     for(j in nrow(p2)){
# 
#       x.sq   <- (p1[i, 1] - p2[j, 1])^2
#       y.sq   <- (p1[i, 2] - p2[j, 2])^2
#       dist.i <- sqrt(x.sq + y.sq)
# 
#       ## test if point is within dist and a different year
#       ## if so, add year to list and check if minimum years has been reached
#       ##        if so break and go to next point
#       if(dist.i <= dist && p2[j,3] != year_i){
# 
#           years <- c(years, p2[j,3])
#           if(length(unique(years)) >= min_years){
#             break
#           }
#       }
#     }
# 
#     ## if there weren't min years nearby, add index to ind_out
#     if(length(unique(years)) < min_years){
#       ind_out <- c(ind_out, i)
#     }
#   }

  return(ind_out)
}



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



auc_pa_gams <- function(gam_list, prediction_data = data_next, target = 'present_total'){
  
  auc <- double(length(gam_list))
  for(i in 1:length(gam_list)){
    
    # Predict presence under future year
    pred_i <- predict(gam_list[[j]], prediction_data, type = 'response', exclude = "s(year_fac)")
    
    # Find AUC metric
    roc_curve <- roc(prediction_data[[target]], pred_j)
    auc_i <- auc(roc_curve)
    
    auc[j] <- auc_i
  }
}



# Validation Functions
#################################################################################

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



# Model Output Functions
#################################################################################

sim_delta_ntimes <- function(gam_bmass, gam_pa, pred_data, n = 5000){
  ###############################
  ## Makes n simulations of delta output given the biomass and presence gams provided.
  ## Returns array with columns containing full prediction sets for a given simulation and 
  ## rows containing all simulated values for a given location
  ## Variables:
  ## gam_bmass is the biomass gam to be used in the simulation
  ## gam_pa is the presence gam to be used in the simulation
  ## pred_data is the prediction dataset used to make the model output
  ## n is the number of output simulations to be created in this function
  
  # Gaussian model simulation
  ###########################################################################################
  # Make predictions on pred_data set
  # Simulate n gaussian output prediction sets using the fit and se.fit
  # - simulating n sets of predictions appends each set in sequence to sim_gauss - tested
  # restructure into array - rows are each point, columns are sets of predictions - tested
  # transform back to normal space using Eric Ward's suggested bias correction
  pred_gauss_with_se <- predict(gam_bmass, newdata = pred_data, type = "link", exclude = "s(year_fac)", se.fit = T)
  pred_gauss_draws <- rnorm(length(pred_gauss_with_se$fit) * n, pred_gauss_with_se$fit, pred_gauss_with_se$se.fit) # randomly selects one value from prob dist of each prediction point based on the given fit (mean) and se
  sim_gauss_array <- array(pred_gauss_draws, dim = c(nrow(pred_data), n))
  sim_gauss <- exp(sim_gauss_array + c((pred_gauss_with_se$se.fit^2)/2)) # bias correction
  
  # Binomial model simulation
  ###########################################################################################
  # Make predictions on pred_data set
  # Simulate n binomial output prediction sets using fit and se.fit
  # restructure just like above section
  # transform binomial to probability space
  pred_binomial_with_se <- predict(gam_pa, newdata = pred_data, type = "link", exclude = "s(year_fac)", se.fit = T)
  pred_binom_draws <- rnorm(length(pred_binomial_with_se$fit) * n, pred_binomial_with_se$fit, pred_binomial_with_se$se.fit)
  sim_binom_array <- array(pred_binom_draws, dim = c(nrow(pred_data), n))
  sim_binomial <- plogis(sim_binom_array)
  
  # Combine predictions
  combined_sim_pred <- sim_gauss * sim_binomial
  
  return(combined_sim_pred)
}



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
  ##  - out_spat_expl - spatially explicit predictions for each prediction year
  
  # Set temporal boundaries for training and prediction in LFO analysis
  start_train_ind <- which(unique(pred_data$year) == start_train)
  years <- unique(pred_data$year)[start_train_ind:length(unique(pred_data$year))]
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
    # over each row, k (each prediction location), find SD and calc CV (sd/mean)
    #################################################################################
    sd_i <- double(length = nrow(simulations_i))
    cv_i <- double(length = nrow(simulations_i))
    for(k in 1:nrow(simulations_i)){
      sd_i[k] <- sd(simulations_i[k, ])
      cv_i[k] <- sd_i[k] / mean_pred_i[k]
    }
    out_spat_expl_list[[i]] <- data.frame(mean = mean_pred_i, 
                                          sd = sd_i,
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
    
    # 1. 
    US_array <- simulations_i[in_US,] 
    CAN_array <- simulations_i[!in_US,]
    
    # 2.
    US_sums <- colSums(US_array)
    CAN_sums <- colSums(CAN_array)
    
    # 3.
    mean_reg <- c(mean(US_sums), mean(CAN_sums))
    sd_reg <- c(sd(US_sums), sd(CAN_sums))
    se_reg <- sd_reg/sqrt(n_sims)
    cv_reg <- sd_reg/mean_reg
    region <- c('US', 'CAN')
    out_reg_list[[i]] <- data.frame(mean_reg, sd_reg, se_reg, cv_reg, region, year = years[i])
  }
  
  out_spat_expl <- do.call(rbind, out_spat_expl_list)
  out_reg <- do.call(rbind, out_reg_list)
  out_list <- list(out_spat_expl, out_reg)
  names(out_list) <- c("out_spat_expl", "out_reg")
  
  return(out_list)
}



# plot functions
###########################################################

plot_heatmap <- function(data, years, age = 'total', names, 
                         type = 'RMSE', lims = c(0,4), data_type = 'J-SCOPE', 
                         xlab = "year predicted"){
  ############################
  # plots accuracy of each tested model in model_selection.R for each year investigated
  # RMSE is plotted for density models
  # AUC is plotted for occurrence models
  
  
  # flip rows for plotting
  dat <- data[rev(1:nrow(data)),]
  
  # Create a heatmap using ggplot2
  dat_df <- melt(dat)
  dat_df$model <- factor(dat_df$Var1)
  dat_df$year <- factor(dat_df$Var2)
  
  if(type == 'RMSE'){
    g <- ggplot(data = dat_df, aes(x = year, y = model, fill = value)) +
      geom_tile() +
      geom_text(aes(label = round(value, 3)), size = 3.5) +
      scale_fill_gradientn(colours = rev(heat.colors(20)),  # You can customize the color palette
                           name = expression(paste('RMSE of\nlog'[e], '(B)')), 
                           limits = lims, 
                           oob = scales::squish) +
      scale_x_discrete(labels = years) +  # Customize X-axis labels
      scale_y_discrete(labels = names)  +    
      theme(axis.text.x=element_text(size=10),
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
      theme(axis.text.x=element_text(size=10),
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
      geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
      geom_point(aes(x = easting_m, y = northing_m, color = mean), shape = 15, size = 0.06) +
      geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
      labs(x = "Longitude",
           y = "Latitude",
           color = 'kg',
           title = paste0('Predicted Biomass for ', years[i])) +
      scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
      ylim(ymin, ymax) +
      scale_color_viridis(limits = delta_lims, oob = scales::squish) +
      theme_simple(base_size = 8, grid = TRUE) +
      theme(plot.title = element_text(hjust = 0.5),
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




plot_ageclass_maps <- function(data, delta_lims, cv_lims, change_lims, age){
  
  #############################
  ## Plots mean forecasted biomass maps, corresponding CV, and change in biomass from 2023 to 2024
  ## data should be of form spatExpl_sim_out--created in generate_output.R
  
  
  plot_delta_pred <- ggplot(data[data$age == age, ]) +
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    geom_point(aes(x = easting_m, y = northing_m, color = sim_delta_2024/1000), shape = 15, size = 0.06) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0('Predicted Biomass')) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_viridis(limits = delta_lims, oob = scales::squish) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(plot.title = element_text(hjust = 0.5),
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
  
  plot_cv <- ggplot(data[data$age == age, ]) +
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    geom_point(aes(x = easting_m, y = northing_m, color = sim_cv_2024), shape = 15, size = 0.06) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
    labs(x = "Longitude",
         y = "Latitude",
         color = '',
         title = paste0('Prediction CV')) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_viridis(limits = cv_lims, oob = scales::squish) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(plot.title = element_text(hjust = 0.5),
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
  # 
  # plot_delta_hist <- ggplot(data[data$age == age,]) +
  #   geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", linewidth = 0.2) +
  #   geom_point(aes(x = easting_m, y = northing_m, color = sim_delta_2023/1000), shape = 15, size = 0.06) +
  #   geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
  #   labs(x = "Longitude",
  #        y = "Latitude",
  #        color = 'kg x 10^3',
  #        title = paste0('Historical Predicted Biomass')) +
  #   scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  #   ylim(ymin, ymax) +
  #   scale_color_viridis(limits = delta_lims, oob = scales::squish) +
  #   theme_simple(base_size = 8, grid = TRUE) +
  #   theme(plot.title = element_text(hjust = 0.5),
  #         panel.spacing.x = unit(-0.5, "pt"),
  #         panel.spacing = unit(1, "lines"),
  #         legend.justification = c(0, 0),
  #         legend.position = c(0.01, 0.005),
  #         legend.title.align = 0.5,
  #         legend.key.size = unit(0.4, "cm"),
  #         legend.background = element_blank(),
  #         legend.text = element_text(size = 6),
  #         legend.title = element_text(size = 6),
  #         panel.grid.major = element_line(size = 0.2))
  # 
  plot_delta_change <- ggplot(data[data$age == age, ]) +
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    geom_point(aes(x = easting_m, y = northing_m, color = delta_change/1000), shape = 15, size = 0.06) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "black", size = 0.2) +  
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0('Change in Biomass')) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_gradient2(low = "steelblue",
                          mid = "grey90",
                          high = "red3",
                          space = "Lab",
                          na.value = "transparent",
                          breaks = c(min(change_lims), min(change_lims)/2, 0 , max(change_lims)/2, max(change_lims)),
                          limits = change_lims,
                          midpoint = 0,
                          oob = scales::squish) +    
    theme_simple(base_size = 8, grid = TRUE) +
    theme(plot.title = element_text(hjust = 0.5),
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
  
  full_plot <- gridExtra::grid.arrange(plot_delta_pred, plot_cv, plot_delta_change, ncol = 3)
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
    theme(plot.background = element_rect(fill = bg, colour = NA),
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


# utility functions
###########################################################


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


rmse <- function(true_vec, predicted_vec){
  ############################
  # find RMSE of predicted vector as compared to the true vector
  if(length(true_vec) != length(predicted_vec)){print('error: trained and predicted vectors not same length')}
  
  return(sqrt(mean((true_vec - predicted_vec)^2, na.rm = T)))
}








