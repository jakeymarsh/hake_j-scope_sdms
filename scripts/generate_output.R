#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Retrospective Reforecasting and 2024 Biomass Prediction
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Generates biomass forecasts using a two-stage Delta Model.
#'              This approach combines binomial (occurrence) and Gaussian 
#'              (density) GAMs to predict spatially explicit biomass signals.
#'              Includes 5,000 Monte Carlo simulations per age class. Compares predicted 
#'              biomass shifts between US and Canadian EEZs.
#' #############################################################################


# ----------------------------------------------------------------------------
# 1. Setup and Data Loading
# ----------------------------------------------------------------------------
# Load selected best model formulas and environmental data
load('./data/temp_js.RData')
load('./data/hake_js.RData')
load('./data/border_chull.RData')
load('./gam/best.formulas.RData')

if(!dir.exists('./outputs/output_stats/')) dir.create('./outputs/output_stats/', recursive = TRUE)

# Global Simulation Parameters
years <- unique(hake_js$year)
nyears <- length(years)
start <- 2012 
num_sims <-  5000
pred_data <- temp_js
train_data <- hake_js

# Create Spatial Grid
grid <- pred_data[pred_data$year == 2015, c("lat", "lon", "easting_m", "northing_m", "read_order")]


# ----------------------------------------------------------------------------
# 2. Run Monte Carlo Delta Model Simulations
# ----------------------------------------------------------------------------

# Total Biomass Forecast
total_sim_region <- lfo_monte_carlo_sim(
  train_data = train_data, 
  pred_data = pred_data, 
  subset = 'wgt_total', 
  start_train = start,
  formula.bmass = best.formulas$bmass.total, 
  formula.pa = best.formulas$pa.total, 
  n_sims = num_sims
)

# Age 2 Biomass Forecast
age2_sim_region <- lfo_monte_carlo_sim(
  train_data = train_data, 
  pred_data = pred_data, 
  subset = 'age2', 
  start_train = start,
  formula.bmass = best.formulas$bmass.age2, 
  formula.pa = best.formulas$pa.age2, 
  n_sims = num_sims
)

# Age 3-4 Biomass Forecast
age3_4_sim_region <- lfo_monte_carlo_sim(
  train_data = train_data, 
  pred_data = pred_data, 
  subset = 'age3_4', 
  start_train = start,
  formula.bmass = best.formulas$bmass.age3_4, 
  formula.pa = best.formulas$pa.age3_4, 
  n_sims = num_sims
)

# Age 5+ Biomass Forecast
age5_20_sim_region <- lfo_monte_carlo_sim(
  train_data = train_data, 
  pred_data = pred_data, 
  subset = 'age5_20', 
  start_train = start,
  formula.bmass = best.formulas$bmass.age5_20, 
  formula.pa = best.formulas$pa.age5_20, 
  n_sims = num_sims
)

save(total_sim_region, file = './outputs/output_stats/total_sim_region.RData')
save(age2_sim_region, file = './outputs/output_stats/age2_sim_region.RData')
save(age3_4_sim_region, file = './outputs/output_stats/age3_4_sim_region.RData')
save(age5_20_sim_region, file = './outputs/output_stats/age5_20_sim_region.RData')



# ----------------------------------------------------------------------------
# 3. Spatially Explicit Simulation Output Analysis
# ----------------------------------------------------------------------------

# Extract annual simulation outputs
spat_expl_simulation_outputs_list <- list()
spat_expl_simulation_outputs_list[[1]] <- data.frame(total_sim_region[["out_spat_expl"]], age = 'Total')
spat_expl_simulation_outputs_list[[2]] <- data.frame(age2_sim_region[["out_spat_expl"]], age = 'Age 2')
spat_expl_simulation_outputs_list[[3]] <- data.frame(age3_4_sim_region[["out_spat_expl"]], age = 'Age 3-4')
spat_expl_simulation_outputs_list[[4]] <- data.frame(age5_20_sim_region[["out_spat_expl"]], age = 'Age 5+')
spat_expl_simulation_outputs <- do.call(rbind, spat_expl_simulation_outputs_list)
spat_expl_simulation_outputs$age <- factor(spat_expl_simulation_outputs$age, levels = c("Total", "Age 2", "Age 3-4", "Age 5+"))



spat_expl_simulation_outputs <- spat_expl_simulation_outputs %>%
  left_join(grid, by = "read_order") %>%
  arrange(age, year, match(read_order, grid$read_order))


save(spat_expl_simulation_outputs, file = './outputs/output_stats/spat_expl_simulation_outputs.RData')



# Make map of biomass change from pre_year to post_year including noting regions of significant changes based on pooled SE

make_change_map <- function(grid, spat_expl_simulation_outputs, pre_year, post_year){
 
  # Calculate mean PRE delta-model forecast map by location/age class using the original row index (read order)
  sim_out_pre_list <- list()
  ages <- unique(spat_expl_simulation_outputs$age)
  nages <- length(ages)
  for(i in 1:nages){
    sim_age_i <- spat_expl_simulation_outputs[spat_expl_simulation_outputs$age == ages[i] & spat_expl_simulation_outputs$year == pre_year,]
    stopifnot(all(grid$read_order == sim_age_i$read_order))  # make sure read_orders line up in grid and simulation
    sim_out_pre_list[[i]] <- data.frame(grid, age = ages[i], year = sim_age_i$year, mean = sim_age_i$mean, cv = sim_age_i$cv, se = sim_age_i$se)
    names(sim_out_pre_list)[i] <- ages[i]
  }
  sim_out_pre <- do.call(rbind, sim_out_pre_list)
  
  # Calculate mean POST delta-model forecast map by location/age class using the original row index (read order)
  sim_out_post_list <- list()
  for(i in 1:nages){
    sim_age_i <- spat_expl_simulation_outputs[spat_expl_simulation_outputs$age == ages[i] & spat_expl_simulation_outputs$year == post_year,]
    stopifnot(all(grid$read_order == sim_age_i$read_order)) # make sure read_orders line up in grid and simulation
    sim_out_post_list[[i]] <- data.frame(grid, age = ages[i], year = sim_age_i$year, mean = sim_age_i$mean, cv = sim_age_i$cv, se = sim_age_i$se)
    names(sim_out_post_list)[i] <- ages[i]
  }
  sim_out_post <- do.call(rbind, sim_out_post_list)
  
  # Ensure spatial alignment between historic and future simulations
  stopifnot(all(sim_out_post$read_order == sim_out_pre$read_order))
  
  # Make final data.frame with locations, ages, historical mean, post forecast
  spatExpl_sim_out <- data.frame(lat = sim_out_pre$lat, 
                                 lon = sim_out_pre$lon, 
                                 easting_m = sim_out_pre$easting_m, 
                                 northing_m = sim_out_pre$northing_m, 
                                 read_order = sim_out_pre$read_order, 
                                 mean_pre = sim_out_pre$mean, 
                                 mean_post = sim_out_post$mean, 
                                 se_pre = sim_out_pre$se, 
                                 se_post = sim_out_post$se, 
                                 cv_pre = sim_out_pre$cv,
                                 cv_post = sim_out_post$cv, 
                                 age = sim_out_post$age)
  spatExpl_sim_out$mean_diff <- spatExpl_sim_out$mean_post - spatExpl_sim_out$mean_pre
  
  # Calculate Pooled SE: sqrt(var1 + var2)
  # Note: Models used to predict 2024 and 2023 share historical training data, making Cov(Y_pre, Y_post) > 0.
  # Assuming independence omits the positive covariance, resulting in a conservative 
  # estimate of variance for the difference (lower Type I error).
  spatExpl_sim_out$pooled_se <- sqrt(spatExpl_sim_out$se_post^2 + spatExpl_sim_out$se_pre^2)
  
  # Calculate Z-score using Pooled SE
  spatExpl_sim_out$change_zscore <- ifelse(
    spatExpl_sim_out$pooled_se > 0,
    spatExpl_sim_out$mean_diff / spatExpl_sim_out$pooled_se,
    0
  )
  
  # Create a mask for non-significant changes for use in plotting (Alpha = 0.05, Z crit = 1.96) 
  spatExpl_sim_out$mean_is_significant <- abs(spatExpl_sim_out$change_zscore) >= 1.96
 
  return(spatExpl_sim_out)
}


change_2023_2024 <- make_change_map(grid, spat_expl_simulation_outputs, 2023, 2024)

spatExpl_sim_out <- change_2023_2024

save(spatExpl_sim_out, file = "./outputs/output_stats/spatExpl_sim_out.RData")



# ----------------------------------------------------------------------------
# 4. Regional Simulation Output Analysis
# ----------------------------------------------------------------------------

# Extract annual regionally aggregated simulation outputs
regional_simulation_outputs_list <- list()
regional_simulation_outputs_list[[1]] <- data.frame(total_sim_region[["out_reg"]], age = 'Total')
regional_simulation_outputs_list[[2]] <- data.frame(age2_sim_region[["out_reg"]], age = 'Age 2')
regional_simulation_outputs_list[[3]] <- data.frame(age3_4_sim_region[["out_reg"]], age = 'Age 3-4')
regional_simulation_outputs_list[[4]] <- data.frame(age5_20_sim_region[["out_reg"]], age = 'Age 5+')
regional_simulation_outputs <- do.call(rbind, regional_simulation_outputs_list)
regional_sim_out <- bind_rows(regional_simulation_outputs)

# Find 95% confidence intervals
regional_sim_out$lower_reg <- pmax(0, regional_sim_out$mean_reg - 1.96 * regional_sim_out$se_reg)
regional_sim_out$upper_reg <- regional_sim_out$mean_reg + 1.96 * regional_sim_out$se_reg
regional_sim_out$age <- factor(regional_sim_out$age, levels = c( "Total", "Age 2", "Age 3-4", "Age 5+"))
save(regional_sim_out, file = './outputs/output_stats/regional_sim_out.RData')


# ----------------------------------------------------------------------------
# 5. Regional Survey Data Analysis
# ----------------------------------------------------------------------------

# Calculate observed hake biomass in US and Canada
ages <- factor(c("Total", "Age 2", "Age 3-4", "Age 5+"), ordered = T)
hake_wgt_obs <- hake_js[ , c("lat", "lon", "year", "wgt_total", "age2", "age3_4", "age5_20")]
hake_obs_reg_list <- list()
for(i in 1:nyears){
  
  # Index for hake in current year
  hake_i <- hake_wgt_obs[hake_wgt_obs$year == years[i],]
  
  # Select hake points in US and CAN
  in_US <- in.chull(hake_i$lon, hake_i$lat, border_chull$lon, border_chull$lat)
  
  # Sum age group biomass in each nation's waters
  sums_US <- colSums(hake_i[in_US, c("wgt_total", "age2", "age3_4", "age5_20")])
  sums_CAN <- colSums(hake_i[!in_US, c("wgt_total", "age2", "age3_4", "age5_20")])
  
  if(i == nyears){
    sums_US <- rep(NA, 4)
    sums_CAN <- rep(NA, 4)
  }
  
  hake_obs_reg_list[[i]] <- data.frame(region = c(rep("US", 4), rep("CAN", 4)),
                                       age = rep(ages, 2),
                                       year = rep(years[i], 8),
                                       wgts = c(sums_US, sums_CAN))
}
hake_obs_reg <- bind_rows(hake_obs_reg_list)
save(hake_obs_reg, file = "./outputs/output_stats/hake_obs_reg.RData")


# end of script


