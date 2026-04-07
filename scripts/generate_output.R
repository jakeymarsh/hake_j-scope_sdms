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
start <- 2012 
num_sims <-  5000
pred_data <- temp_js
train_data <- hake_js


# ----------------------------------------------------------------------------
# 2. Run Monte Carlo Delta Model Simulations
# ----------------------------------------------------------------------------

# Total Biomass Forecast
total_sim_region <- lfo_monte_carlo_sim(
  train_data = train_data, 
  pred_data = pred_data, 
  subset = 'wgt_total', 
  start_train = start,
  formula.bmass = best.formulas$bmass.total, # Stage 2: Density
  formula.pa = best.formulas$pa.total, # Stage 1: Occurrence
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
save(spat_expl_simulation_outputs, file = './outputs/output_stats/spat_expl_simulation_outputs.RData')



# Make calculate map of biomass change from pre_year to post_year including noting regions of significant changes based on pooled SE

make_change_map <- function(grid, spat_expl_simulation_outputs, pre_year, post_year){
 
  # Calculate mean pre delta-model forecast map by location/age class using the original row index (read order)
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
  
  # Calculate post mean delta-model forecast map by location/age class using the original row index (read order)
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
  
  # Calculate Pooled SE: sqrt(var1 + var2) -- (assumes forecasts of the two years are independent)
  spatExpl_sim_out$pooled_se <- sqrt(spatExpl_sim_out$se_post^2 + spatExpl_sim_out$se_pre^2)
  
  # Calculate Z-score using Pooled SE
  spatExpl_sim_out$change_zscore <- spatExpl_sim_out$mean_diff / spatExpl_sim_out$pooled_se
  
  # Create a mask for non-significant changes for use in plotting (Alpha = 0.05, Z crit = 1.96) 
  spatExpl_sim_out$mean_is_significant <- abs(spatExpl_sim_out$change_zscore) >= 1.96
 
  return(spatExpl_sim_out)
}

grid <- pred_data[pred_data$year == 2015, c("lat", "lon", "easting_m", "northing_m", "read_order")]

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
regional_sim_out$lower <- regional_sim_out$mean - 1.96 * regional_sim_out$se
regional_sim_out$upper <- regional_sim_out$mean + 1.96 * regional_sim_out$se
regional_sim_out$age <- factor(regional_sim_out$age, levels = c( "Total", "Age 2", "Age 3-4", "Age 5+"))
save(regional_sim_out, file = './outputs/output_stats/regional_sim_out.RData')


# Find predicted diff between US and CAN hake
years <- unique(regional_sim_out$year)
nyears <- length(years)
ages <- unique(regional_sim_out$age)
count <- 1
diffs <- list()
for(i in 1:length(years)){
  for(j in 1:length(ages)){
    means_i <- regional_sim_out$mean[regional_sim_out$year == years[i] & regional_sim_out$age == ages[j]]
    diffs[[count]] <- data.frame(diff = -diff(means_i),
                                 sum = sum(means_i),
                                 ratio_US = -diff(means_i) / sum(means_i),
                                 year = years[i],
                                 age = ages[j])
    count <- count+1
  }
}
region_diff <- bind_rows(diffs)
save(region_diff, file = './outputs/output_stats/region_diff.RData')


# ----------------------------------------------------------------------------
# 5. Regional Survey Data Analysis
# ----------------------------------------------------------------------------

# Calculate total biomass in US and CAN, as observed on hake surveys
years <- unique(hake_obs_reg$year)
ages <- unique(hake_obs_reg$age)
count <- 1
hake_obs_reg_diff_list <- list()
for(i in 1:length(years)){
  for(j in 1:length(ages)){
    wgts_i <- hake_obs_reg$wgts[hake_obs_reg$year == years[i] & hake_obs_reg$age == ages[j]]
    hake_obs_reg_diff_list[[count]] <- data.frame(diff = -diff(wgts_i),
                                          sum = sum(wgts_i),
                                          ratio_US = -diff(wgts_i) / sum(wgts_i),
                                          year = years[i],
                                          age = ages[j])
    count <- count+1
  }
}
hake_obs_reg_diff_list[[count]] <- data.frame(diff = rep(NA, 4),
                                      sum = rep(NA, 4),
                                      ratio_US = rep(NA, 4),
                                      year = 2024,
                                      age = ages)
hake_obs_reg_diff <- bind_rows(hake_obs_reg_diff_list)
save(hake_obs_reg_diff, file = './outputs/output_stats/hake_obs_reg_diff.RData')


# end of script


