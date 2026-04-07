######################################################################
## Generate model forecasts and reforecasts for retrospective time period (2015 - 2023)
## Forecast spatially explicitly (maps) and regionally aggregated outputs
######################################################################



######################################################################
## Use lfo_monte_carlo_sim() to predict biomass over retrospective period and 2024
## - using best models for each age class, selected in model_selection.R
## - using 5000 simulations, only looking at 2015 onward
######################################################################
# load hake data, J-SCOPE forecasts, and border convex hull
load('./data/temp_js.RData')
load('./data/hake_js.RData')
load('./data/border_chull.RData')
load('./gam/best.formulas.RData')


if(!dir.exists('./outputs/output_stats/'))
  {dir.create('./outputs/output_stats/')}

start <- 2012 
num_sims <-  5000
pred_data <- temp_js
train_data <- hake_js

total_sim_region <- lfo_monte_carlo_sim(train_data = train_data, 
                                                   pred_data = pred_data, 
                                                   subset = 'wgt_total', 
                                                   start_train = start,
                                                   formula.bmass = best.formulas$bmass.total, 
                                                   formula.pa = best.formulas$pa.total, 
                                                   n_sims = num_sims)

age2_sim_region <- lfo_monte_carlo_sim(train_data = train_data, 
                                                  pred_data = pred_data, 
                                                  subset = 'age2', 
                                                  start_train = start,
                                                  formula.bmass = best.formulas$bmass.age2, 
                                                  formula.pa = best.formulas$pa.age2, 
                                                  n_sims = num_sims)

age3_4_sim_region <- lfo_monte_carlo_sim(train_data = train_data, 
                                                    pred_data = pred_data, 
                                                    subset = 'age3_4', 
                                                    start_train = start,
                                                    formula.bmass = best.formulas$bmass.age3_4, 
                                                    formula.pa = best.formulas$pa.age3_4, 
                                                    n_sims = num_sims)

age5_20_sim_region <- lfo_monte_carlo_sim(train_data = train_data, 
                                                     pred_data = pred_data, 
                                                     subset = 'age5_20', 
                                                     start_train = start,
                                                     formula.bmass = best.formulas$bmass.age5_20, 
                                                     formula.pa = best.formulas$pa.age5_20, 
                                                     n_sims = num_sims)


###########################################################################
## Spatially Explicit Simulation Output Analysis
###########################################################################

# extract annual simulation outputs
spat_expl_simulation_outputs_list <- list()
spat_expl_simulation_outputs_list[[1]] <- data.frame(total_sim_region[["out_spat_expl"]], age = 'Total')
spat_expl_simulation_outputs_list[[2]] <- data.frame(age2_sim_region[["out_spat_expl"]], age = 'Age 2')
spat_expl_simulation_outputs_list[[3]] <- data.frame(age3_4_sim_region[["out_spat_expl"]], age = 'Age 3-4')
spat_expl_simulation_outputs_list[[4]] <- data.frame(age5_20_sim_region[["out_spat_expl"]], age = 'Age 5+')
spat_expl_simulation_outputs <- do.call(rbind, spat_expl_simulation_outputs_list)
spat_expl_simulation_outputs$age <- factor(spat_expl_simulation_outputs$age, levels = c("Total", "Age 2", "Age 3-4", "Age 5+"))
save(spat_expl_simulation_outputs, file = './outputs/output_stats/spat_expl_simulation_outputs.RData')

# finds mean of all points with same read_order (same location) taken from grid formed here
grid <- pred_data[pred_data$year == 2015, c("lat", "lon", "easting_m", "northing_m", "read_order")]

# find mean delta output of 2023 for each location using with read_order (location-specific index) for each age class
hist_sim_list <- list()
ages <- unique(spat_expl_simulation_outputs$age)
nages <- length(ages)
for(i in 1:nages){
  sim_age_i <- spat_expl_simulation_outputs[spat_expl_simulation_outputs$age == ages[i] & spat_expl_simulation_outputs$year == 2023,]
  hist_sim_i <- tapply(sim_age_i$mean, as.factor(sim_age_i$read_order), mean)
  hist_sim_list[[i]] <- data.frame(grid, age = ages[i], sim_delta_hist = hist_sim_i)
  names(hist_sim_list)[i] <- ages[i]
}
hist_sim <- do.call(rbind, hist_sim_list)

# find 2024 delta forecast for each age class
fut_sim_list <- list()
for(i in 1:nages){
  sim_age_i <- spat_expl_simulation_outputs[spat_expl_simulation_outputs$age == ages[i] & spat_expl_simulation_outputs$year == 2024,]
  fut_sim_list[[i]] <- data.frame(grid, age = ages[i], mean = sim_age_i$mean, cv = sim_age_i$cv)
  names(fut_sim_list)[i] <- ages[i]
}
fut_sim <- do.call(rbind, fut_sim_list)

# make sure read orders match up
if(any(hist_sim$read_order != fut_sim$read_order)){print('means not matching up spatially, fix this!')}

# make final data.frame with locations, ages, historical mean, 2024 forecast
spatExpl_sim_out <- data.frame(lat = hist_sim$lat, 
                               lon = hist_sim$lon, 
                               easting_m = hist_sim$easting_m, 
                               northing_m = hist_sim$northing_m, 
                               read_order = hist_sim$read_order, 
                               sim_delta_2023 = hist_sim$sim_delta_hist, 
                               sim_delta_2024 = fut_sim$mean, 
                               sim_cv_2024 = fut_sim$cv, 
                               age = hist_sim$age)
spatExpl_sim_out$delta_change <- spatExpl_sim_out$sim_delta_2024 - spatExpl_sim_out$sim_delta_2023

# take out all changes that are not significantly larger than 0
zscore <- spatExpl_sim_out$delta_change/spatExpl_sim_out$sim_se_2024
spatExpl_sim_out$delta_change[abs(zscore)<1.96] <- NA
spatExpl_sim_out$percent_change <- spatExpl_sim_out$delta_change/spatExpl_sim_out$sim_delta_2023*100

save(spatExpl_sim_out, file = "./outputs/output_stats/spatExpl_sim_out.RData")


###########################################################################
## Regional Simulation Output Analysis
###########################################################################

# extract annual regionally aggregated simulation outputs
regional_simulation_outputs_list <- list()
regional_simulation_outputs_list[[1]] <- data.frame(total_sim_region[["out_reg"]], age = 'Total')
regional_simulation_outputs_list[[2]] <- data.frame(age2_sim_region[["out_reg"]], age = 'Age 2')
regional_simulation_outputs_list[[3]] <- data.frame(age3_4_sim_region[["out_reg"]], age = 'Age 3-4')
regional_simulation_outputs_list[[4]] <- data.frame(age5_20_sim_region[["out_reg"]], age = 'Age 5+')
regional_simulation_outputs <- do.call(rbind, regional_simulation_outputs_list)
regional_sim_out <- bind_rows(regional_simulation_outputs)

# find 95% confidence intervals
regional_sim_out$lower <- regional_sim_out$mean - 1.96 * regional_sim_out$se
regional_sim_out$upper <- regional_sim_out$mean + 1.96 * regional_sim_out$se
regional_sim_out$age <- factor(regional_sim_out$age, levels = c( "Total", "Age 2", "Age 3-4", "Age 5+"))
save(regional_sim_out, file = './outputs/output_stats/regional_sim_out.RData')


# find predicted diff between US and CAN hake
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





##########################################################
# Regional Survey Data Analysis
##########################################################

# Ratio of hake in US
ages <- factor(c("Total", "Age 2", "Age 3-4", "Age 5+"), ordered = T)
hake_wgt_obs <- hake_js[ , c("lat", "lon", "year", "wgt_total", "age2", "age3_4", "age5_20")]
hake_obs_reg_list <- list()
for(i in 1:nyears){
  
  # find hake in current year
  hake_i <- hake_wgt_obs[hake_wgt_obs$year == years[i],]
  
  # select hake points in US and CAN
  in_US <- in.chull(hake_i$lon, hake_i$lat, border_chull$lon, border_chull$lat)
  
  # sum them up for each age class and find ratio in the US waters
  sums_US <- colSums(hake_i[in_US, c("wgt_total", "age2", "age3_4", "age5_20")])
  sums_CAN <- colSums(hake_i[!in_US, c("wgt_total", "age2", "age3_4", "age5_20")])
  ratio_US <- sums_US / colSums(matrix(c(sums_US, sums_CAN), nrow = 2, ncol = 4, byrow = T))
  
  # find ratio of # observations in US to total
  num_obs_US <- length(which(in_US))
  num_obs_CAN <- length(which(!in_US))
  ratio_num_obs_US <- num_obs_US / (num_obs_US + num_obs_CAN)
  
  
  if(i == nyears){
    sums_US <- rep(NA, 4)
    sums_CAN <- rep(NA, 4)
  }
  
  hake_obs_reg_list[[i]] <- data.frame(region = c(rep("US", 4), rep("CAN", 4)),
                                       age = rep(ages, 2),
                                       year = rep(years[i], 8),
                                       wgts = c(sums_US, sums_CAN), 
                                       ratio_US = rep(ratio_US, 2),
                                       ratio_num_obs_US = ratio_num_obs_US)
}
hake_obs_reg <- bind_rows(hake_obs_reg_list)
save(hake_obs_reg, file = "./outputs/output_stats/hake_obs_reg.RData")


# Biomass in US and CAN
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




##########################################################
# Inter Annual Temperature Forecast Analysis
##########################################################

# find 2024 predicted temperature
pred_data <- temp_js
pred_2024<- pred_data[pred_data$year == 2024,]
hist_fut_temp <- data.frame(t250_fut = pred_2024$t250m_Ano_fill,
                            ind = 1:nrow(pred_2024),
                            easting_m = pred_2024$easting_m, 
                            northing_m = pred_2024$northing_m,
                            lat=  pred_2024$lat, 
                            lon =  pred_2024$lon)

# add in 2023 temperature
years <- c(2012, 2013, 2015, 2017, 2019, 2021, 2023)
nyears <- length(years)
hist_temp <- pred_data[pred_data$year == 2023, "t250m_Ano_fill"]

# consolidate data in hist_fut_temp
hist_fut_temp$t250_hist <- hist_temp
hist_fut_temp$t250_change <-  hist_fut_temp$t250_fut - hist_fut_temp$t250_hist
hist_fut_temp$percent_change <-  hist_fut_temp$t250_change / hist_fut_temp$t250_hist * 100

save(hist_fut_temp, file = "./outputs/output_stats/hist_fut_temp.RData")



# end of script


