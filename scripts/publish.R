#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Figures for Publication
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Generates publication-quality figures including study area maps, 
#'              raw survey data visualizations, LFO-CV performance heatmaps, 
#'              partial effects plots, and 2024 biomass forecasts.
#' #############################################################################

# --------------------------------------------------------------------------------
# 1. Workspace Configuration
# --------------------------------------------------------------------------------

# 1.1 Define all required output subdirectories
figure_dirs <- c(
  "./final_outputs/figures/data/",
  "./final_outputs/figures/model_performance/total/",
  "./final_outputs/figures/model_performance/age2/",
  "./final_outputs/figures/model_performance/age3_4/",
  "./final_outputs/figures/model_performance/age5_20/",
  "./final_outputs/figures/model_partial_effects/total/",
  "./final_outputs/figures/model_partial_effects/age2/",
  "./final_outputs/figures/model_partial_effects/age3_4/",
  "./final_outputs/figures/model_partial_effects/age5_20/",
  "./final_outputs/figures/forecasts/age_class_maps/"
)
invisible(lapply(figure_dirs, function(d) if(!dir.exists(d)) dir.create(d, recursive = TRUE)))

# 1.2 Global Plotting Constants (Hyperparameters)
cv_lims    <- c(0, 0.5)  # Prediction Uncertainty (CV)
mean_lims  <- c(0, 3e4)  # Mean Biomass (kg)
delta_lims <- c(0, 5e4)  # Net Biomass Change
rmse_lims  <- c(1, 3)    # LFO-CV Density Error
auc_lims   <- c(0.5, 1)  # LFO-CV Occurrence Accuracy
n_models <- 5

# 1.3 Load Analysis Assets

# 1.4 Model Performance Data
load('./outputs/model_performance/performance_7year_array.RData')
load('./outputs/model_performance/model.performance.RData')

# 1.5 Forecast & Regional Simulation Statistics
load('./outputs/output_stats/spat_expl_simulation_outputs.RData')
load('./outputs/output_stats/spatExpl_sim_out.RData')
load('./outputs/output_stats/regional_sim_out.RData')
load('./outputs/output_stats/region_diff.RData')
load('./outputs/output_stats/hake_obs_reg.RData')
load('./outputs/output_stats/hist_fut_temp.RData')

# 1.6 Best Models, Survey Data, and Oceanographic Forecasts
load('./gam/best.gams.RData')
load('./data/hake_js.RData')
hake_raw_data <- create_vars(read.csv("./data/hake_biomass_2024_01_12.csv"))


#1.7 Define global plotting limits (UTM Zone 8)
lims <- c(0,17)
xmin <- min(hake_js[["easting_m"]], na.rm = TRUE) - 5e3
xmax <- max(hake_js[["easting_m"]], na.rm = TRUE) + 10e3
ymin <- min(hake_js[["northing_m"]], na.rm = TRUE) - 5e3
ymax <- max(hake_js[["northing_m"]], na.rm = TRUE) + 5e3

dat_t <- hake_js[, c("easting_m", "northing_m", "lat", "lon", "shelf_distance_m", "year", "year_fac", "log_wgt_total", "wgt_total")]
dat_2 <- hake_js[, c("easting_m", "northing_m", "lat", "lon", "shelf_distance_m", "year", "year_fac", "log_age2", "age2")]
dat_3 <- hake_js[, c("easting_m", "northing_m", "lat", "lon", "shelf_distance_m", "year", "year_fac", "log_age3_4", "age3_4")]
dat_5 <- hake_js[, c("easting_m", "northing_m", "lat", "lon", "shelf_distance_m", "year", "year_fac", "log_age5_20", "age5_20")]
nam <- c("easting_m", "northing_m", "lat", "lon", "shelf_distance_m", "year", "year_fac", "log_bmass", "bmass")

names(dat_t) <- nam
names(dat_2) <- nam
names(dat_3) <- nam
names(dat_5) <- nam

dat_t$age <- "Total"
dat_2$age <- "Age 2"
dat_3$age <- "Age 3-4"
dat_5$age <- "Age 5+"

hake_dat <- rbind(dat_t, dat_2, dat_3, dat_5)
hake_dat$age <- factor(hake_dat$age, levels = c("Total", "Age 2", "Age 3-4", "Age 5+")) ## data.frame with all hake_js data


# --------------------------------------------------------------------------------
# 2. Map of Region of Study 
# --------------------------------------------------------------------------------

o1 <- ggplot() +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          linewidth = 0.2) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax= ymax), 
            color = "red", fill = NA) +
  scale_x_continuous(breaks = c(seq(-200, -50, 15)), limits = c(xmin-500000*2, xmax+700000*2)) +
  scale_y_continuous(breaks = c(seq(20, 70, 10)), limits = c(ymin-1000000*2, ymax+1000000*1.5)) +
  labs(x = "Longitude",
       y = "Latitude") +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))+
  theme_bw(base_size = 10, base_family = "serif")

ggsave('./final_outputs/figures/placeofstudy.png', o1,
       width = 3, height = 4, units = "in", bg = "white")

ggsave('./final_outputs/figures/placeofstudy.tiff', o1,
       device = "tiff", width = 3, height = 4, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# --------------------------------------------------------------------------------
# 3. Hake Data 
# --------------------------------------------------------------------------------

# Observed biomass over all years for each age class
h1 <- ggplot(hake_dat) +
  geom_point(aes(x = easting_m, y = northing_m, color = log_bmass), 
             shape = 15, size = 0.4, alpha = 0.2) +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          size = 0.2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  scale_color_gradient(low = "#E0E6FA", high = "#3A00BB",  
                       na.value = "transparent", limits = c(0,17)) +
  facet_wrap(~ age, nrow = 1) + 
  scale_x_continuous(breaks = c(-128, -126, -124), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = "Longitude",
       y = "Latitude",
       color = expression(log[e](B))) +
  theme_simple(base_size = 8, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/hake_by_ageclass.png', h1,
       width = 7, height = 3, units = "in", bg = "white")

ggsave('./final_outputs/figures/data/hake_by_ageclass.tiff', h1,
       device = "tiff", width = 7, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# Age 2 hake biomass by year plotted latitude against distance to shelf
h1_2 <- ggplot(hake_raw_data) +
  geom_hline(aes(yintercept = 49.95222), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_hline(aes(yintercept = 42.95613), linetype = "dashed", color = "red", alpha = 0.7) +
  geom_vline(aes(xintercept = 0), linetype = "solid", alpha = 0.7) +
  geom_point(aes(x = shelf_distance_m/1000, y = lat, color = log_age2), size = 0.001) +
  scale_color_viridis(na.value = "transparent") +
  scale_x_reverse() +
  scale_y_continuous(limits = c(34,52))+
  facet_wrap(~year, nrow = 4) +
  labs(x = "Shelf Distance (km)",
       y = "Latitude (ºN)", 
       color = expression(log[e](B)))+
  theme_simple(base_size = 10, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.85, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/age2_lat_shelfdist.png', h1_2,
       width = 6, height = 8, units = "in", bg = "white")

ggsave('./final_outputs/figures/data/age2_lat_shelfdist.tiff', h1_2,
       device = "tiff", width = 6, height = 8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# Hake biomass by year plotted latitude against distance to shelf
h1_3 <- ggplot(hake_js) +
  geom_hline(aes(yintercept = 48.25), linetype = "dashed", alpha = 0.7) +
  geom_vline(aes(xintercept = 0), linetype = "solid", alpha = 0.7) +
  geom_point(aes(x = shelf_distance_m/1000, y = lat, color = log_wgt_total), size = 0.001) +
  scale_color_viridis(na.value = "transparent") +
  scale_x_reverse() +
  facet_wrap(~year_fac, nrow = 4) +
  labs(x = "Shelf Distance (km)",
       y = "Latitude (ºN)", 
       color = expression(log[e](B)))+
  theme_simple(base_size = 10, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.85, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/hake_lat_shelfdist.png', h1_3,
       width = 6, height = 8, units = "in", bg = "white")

ggsave('./final_outputs/figures/data/hake_lat_shelfdist.tiff', h1_3,
       device = "tiff", width = 6, height = 8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --------------------------------------------------------------------------------
# 4. Temperature Forecasts
# --------------------------------------------------------------------------------

# Prep data for plotting
data <- temp_js
pd_t250 <- data[ , c("easting_m", "northing_m", "year_fac", "t250m_Ano_fill")]
pd_t100 <- data[ , c("easting_m", "northing_m", "year_fac", "t100m_Ano_fill")]
nam <- c("easting_m", "northing_m", "year_fac", "temp_anom")
names(pd_t250) <- nam
names(pd_t100) <- nam
pd_t250$temp_depth <- "T250"
pd_t100$temp_depth <- "T100"
pd <- rbind(pd_t250, pd_t100)
pd$temp_depth <- factor(pd$temp_depth, levels = c("T250", "T100"))

# Plot Temperature anomaly at 250m depth for all years
t1 <- ggplot(pd[pd$temp_depth == "T250",]) +
  geom_point(aes(x = easting_m, y = northing_m, color = temp_anom), shape = 15, size = 0.2) +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          size = 0.2) +
  facet_wrap(~year_fac, nrow = 2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1,1),
                        midpoint = 0,
                        oob = scales::squish) +  
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = "Longitude",
       y = "Latitude",
       color = expression(paste("ºC")), 
       title = "Temperature Anomaly at 250 m") +
  theme_simple(base_size = 10, grid = TRUE) +
theme(text = element_text(family = "serif"),
      axis.text = element_text(size = 8, family = "serif"),
      axis.title = element_text(size = 10, family = "serif"),
      legend.text = element_text(size = 8, family = "serif"), 
      legend.title = element_text(size = 10, family = "serif"),
      panel.spacing.x = unit(-0.5, "pt"),
      panel.spacing = unit(1, "lines"),
      legend.justification = c(0, 0),
      legend.position = c(0.9, 0.005),
      legend.title.align = 0.5,
      legend.key.size = unit(0.4, "cm"),
      legend.background = element_blank(),
      panel.grid.major = element_line(size = 0.2),
      strip.text = element_text(family = "serif", size = 10, face = "plain"),
      strip.background = element_blank())

ggsave('./final_outputs/figures/data/js_t250_all_years.png', t1,
       width = 7, height = 4.5, units = "in", bg = "white")
ggsave('./final_outputs/figures/data/js_t250_all_years.tiff', t1,
       device = "tiff", width = 7, height = 4.5, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Plot Temperature anomaly at 100m depth for all years
t2 <- ggplot(pd[pd$temp_depth == "T100",]) +
  geom_point(aes(x = easting_m, y = northing_m, color = temp_anom), shape = 15, size = 0.2) +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          size = 0.2) +
  facet_wrap(~year_fac, nrow = 2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1,1),
                        midpoint = 0,
                        oob = scales::squish) +  
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = "Longitude",
       y = "Latitude",
       color = expression(paste("ºC")), 
       title = "Temperature Anomaly at 100 m") +
  theme_simple(base_size = 10, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.9, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/js_t100_all_years.png', t2,
       width = 7, height = 4.5, units = "in", bg = "white")
ggsave('./final_outputs/figures/data/js_t100_all_years.tiff', t2,
       device = "tiff", width = 7, height = 4.5, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Plot examples of a hot year forecast, a cold year forecast, and the 2024 forecast 
all_years <- unique(temp_js$year)
all_nyears <- length(all_years)
list_mean_temp <- list()
for(i in 1:all_nyears){
  mean_t250_i <- mean(temp_js[temp_js$year == all_years[i], 't250m_Ano_fill'], na.rm = T)
  mean_t100_i <- mean(temp_js[temp_js$year == all_years[i], 't100m_Ano_fill'], na.rm = T)
  list_mean_temp[[i]] <- data.frame(t250 = mean_t250_i, t100 = mean_t100_i, year = all_years[i])
}
mean_temps <- bind_rows(list_mean_temp, .id = "column_label")

# hot and cold are the same for both t100 and t250
cold_ind <- order(mean_temps$t250, decreasing = F)[1] # coldest year, averaged across space
hot_ind <- order(mean_temps$t250, decreasing = F)[nrow(mean_temps)] # hottest year, as averaged across space
fut_ind <- which(mean_temps$year == 2024) 
years_hcf <- mean_temps$year[c(hot_ind, cold_ind, fut_ind)]

dat_hcf <- hake_dat[is.element(hake_dat$year_fac, c(years_hcf)), ]
dat_hcf$year_fac <- factor(dat_hcf$year_fac, levels = c(years_hcf)) # data.frame with only hot, neutral, and cold years


# Plot t250 over hot/cold/2024 years
####################################################################
pd_hcf <- pd[is.element(pd$year_fac, c(years_hcf)), ]
pd_hcf$year_fac <- factor(pd_hcf$year, levels = c(years_hcf))

t4 <- ggplot(pd_hcf[which(pd_hcf$temp_depth == "T250"),]) +
  geom_point(aes(x = easting_m, y = northing_m, color = temp_anom), shape = 15, size = 0.2) +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          size = 0.2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  facet_wrap(~ year_fac, nrow = 1) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1,1),
                        midpoint = 0,
                        oob = scales::squish) +  
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = "Longitude",
       y = "Latitude",
       color = expression(paste("ºC")), 
       title = "Temperature Anomaly at 250 m") +
  theme_simple(base_size = 10, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/js_t250_hot_cold_fut_years.png', t4,
       width = 4, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/data/js_t250_hot_cold_fut_years.tiff', t4,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Plot t100 over hot/cold/2024 years
####################################################################
t5 <- ggplot(pd_hcf[which(pd_hcf$temp_depth == "T100"),]) +
  geom_point(aes(x = easting_m, y = northing_m, color = temp_anom), shape = 15, size = 0.2) +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          size = 0.2) +
  facet_wrap(~ year_fac, nrow = 1) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1,1),
                        midpoint = 0,
                        oob = scales::squish) +  
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = "Longitude",
       y = "Latitude",
       color = expression(paste("ºC")), 
       title = "Temperature Anomaly at 100 m") +
  theme_simple(base_size = 10, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 8, family = "serif"),
        axis.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif"), 
        legend.title = element_text(size = 10, family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        strip.text = element_text(family = "serif", size = 10, face = "plain"),
        strip.background = element_blank())

ggsave('./final_outputs/figures/data/js_t100_hot_cold_fut_years.png', t5,
       width = 4, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/data/js_t100_hot_cold_fut_years.tiff', t5,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# --------------------------------------------------------------------------------
# 5. LFO-CV Model Performance Heatmaps (RMSE and AUC) 
# --------------------------------------------------------------------------------

# Create plotting variables 
n_models <- nrow(performance_7year_array$total_bmass)
model_names <- paste0("Model ", n_models:1)
pred_years_7years <- c(tail(unique(hake_js$year), 7), 'mean')

# Total Density
mp_t_1 <- plot_heatmap(data = performance_7year_array$total_bmass, years = pred_years_7years, 
                       names = model_names, type = 'RMSE', lims = rmse_lims, 
                       age = "Total")
ggsave("./final_outputs/figures/model_performance/total/total_bmass_RMSE_7years.png", mp_t_1,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/total/total_bmass_RMSE_7years.tiff", mp_t_1,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Total Occurrence
mp_t_2 <- plot_heatmap(data = performance_7year_array$total_pa, years = pred_years_7years, 
                       names = model_names, type = 'AUC', lims = auc_lims, 
                       age = "Total")
ggsave("./final_outputs/figures/model_performance/total/total_AUC_7years.png", mp_t_2,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/total/total_AUC_7years.tiff", mp_t_2,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Age 2 Density
mp_a2_1 <- plot_heatmap(data = performance_7year_array$age2_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 2")
ggsave("./final_outputs/figures/model_performance/age2/age2_bmass_RMSE_7years.png", mp_a2_1,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age2/age2_bmass_RMSE_7years.tiff", mp_a2_1,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Age 2 Occurrence
mp_a2_2 <- plot_heatmap(data = performance_7year_array$age2_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 2")
ggsave("./final_outputs/figures/model_performance/age2/age2_AUC_7years.png", mp_a2_2,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age2/age2_AUC_7years.tiff", mp_a2_2,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# Age 3-4 Density
mp_a3_1 <- plot_heatmap(data = performance_7year_array$age3_4_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 3-4")
ggsave("./final_outputs/figures/model_performance/age3_4/age3_4_bmass_RMSE_7years.png", mp_a3_1,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age3_4/age3_4_bmass_RMSE_7years.tiff", mp_a3_1,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Age 3-4 Occurrence
mp_a3_2 <- plot_heatmap(data = performance_7year_array$age3_4_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 3-4")
ggsave("./final_outputs/figures/model_performance/age3_4/age3_4_AUC_7years.png", mp_a3_2,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age3_4/age3_4_AUC_7years.tiff", mp_a3_2,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Age 5+ Density
mp_a5_1 <- plot_heatmap(data = performance_7year_array$age5_20_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 5+")
ggsave("./final_outputs/figures/model_performance/age5_20/age5_20_bmass_RMSE_7years.png", mp_a5_1,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age5_20/age5_20_bmass_RMSE_7years.tiff", mp_a5_1,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# Age 5+ Occurrence
mp_a5_2 <- plot_heatmap(data = performance_7year_array$age5_20_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 5+")
ggsave("./final_outputs/figures/model_performance/age5_20/age5_20_AUC_7years.png", mp_a5_2,
       width = 6, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_performance/age5_20/age5_20_AUC_7years.tiff", mp_a5_2,
       device = "tiff", width = 6, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --------------------------------------------------------------------------------
# 6. Fully Trained Model Performance
# --------------------------------------------------------------------------------

# Standardize the row count across performance lists to ensure clean merge
for(i in 1:length(model.performance)){
  model.performance[[i]] <- model.performance[[i]][1:nrow(model.performance$total.bmass),]
}

# Combine list elements into a single long-format dataframe for display
mp <- bind_rows(model.performance)

# Align occurrence (AUC) and density (RMSE) metrics into a single cohesive table
mp$avg_auc <- mp$avg_auc[!is.na(mp$avg_auc)]
names(mp) <- c("bmass.aic", "bmass.r2", "bmass.dev.expl", "model", "age", "avg_rmse", "avg_auc") #, "pa.aic")

# Formatting: only display 'Age' label on first row of each group for visual clarity
mp$age <- ifelse((1:nrow(mp)) %% n_models == 1, mp$age, NA)
mp <- mp %>% filter(!is.na(avg_rmse))


# Calculate Relative Performance (Delta Metrics) ----------------------------
# We calculate Delta (Δ) values to identify the best-performing model 
# (lowest RMSE or highest AUC) within each specific age-class group

# --- Δ RMSE (Density Precision) ---
# Lower values indicate higher precision. ΔRMSE = 0 represents the 'best' model
delta_rmse <- double(length = nrow(mp))

# Calculate relative error for each of the four hake age classes
delta_rmse[1:n_models] <- mp$avg_rmse[1:n_models] - min(mp$avg_rmse[1:n_models])
delta_rmse[1*n_models+(1:n_models)] <- mp$avg_rmse[1*n_models+(1:n_models)] - min(mp$avg_rmse[1*n_models+(1:n_models)])
delta_rmse[2*n_models+(1:n_models)] <- mp$avg_rmse[2*n_models+(1:n_models)] - min(mp$avg_rmse[2*n_models+(1:n_models)])
delta_rmse[3*n_models+(1:n_models)] <- mp$avg_rmse[3*n_models+(1:n_models)] - min(mp$avg_rmse[3*n_models+(1:n_models)])

# --- Δ AUC (Occurrence Accuracy) ---
# Higher values indicate better discrimination. ΔAUC = 0 represents the 'Best' model
delta_auc <- double(length = nrow(mp))

# Calculate relative accuracy for each of the four hake age classes
delta_auc[1:n_models] <- mp$avg_auc[1:n_models] - max(mp$avg_auc[1:n_models])  
delta_auc[1*n_models+(1:n_models)] <- mp$avg_auc[1*n_models+(1:n_models)] - max(mp$avg_auc[1*n_models+(1:n_models)]) 
delta_auc[2*n_models+(1:n_models)] <-  mp$avg_auc[2*n_models+(1:n_models)] - max(mp$avg_auc[2*n_models+(1:n_models)])
delta_auc[3*n_models+(1:n_models)] <- mp$avg_auc[3*n_models+(1:n_models)] - max(mp$avg_auc[3*n_models+(1:n_models)])

# Prepare the final subset of metrics to be passed to the gt_table for visualization.
mp$delta_rmse <- delta_rmse
mp$delta_auc <- delta_auc
mp_plot <- mp[, c("age", "model", "avg_rmse", "delta_rmse", "bmass.r2", "avg_auc","delta_auc")]

# Create the table
gt_table <- gt(mp_plot) %>%
  tab_header(
    title = "Model Performance"
  ) %>%
  cols_label(
    age = "Age",
    model = "Model",
    avg_rmse = "RMSE",
    delta_rmse = "ΔRMSE",
    bmass.r2 = "R²",
    avg_auc = "AUC",
    delta_auc = "ΔAUC",
  ) %>%
  fmt_number(
    columns = c(avg_rmse, delta_rmse, bmass.r2, avg_auc, delta_auc),
    decimals = 3
  ) %>%
  tab_spanner(
    label = "Density Metrics",
    columns = c(avg_rmse, delta_rmse, bmass.r2),
  ) %>%
  tab_spanner(
    label = "Occurrence Metrics",
    columns = c(avg_auc, delta_auc),
  ) %>%
  tab_options(
    table.width = pct(60),
    column_labels.font.weight = "bold",
    data_row.padding = px(5),
    table.font.size = 16
  ) %>% sub_missing(
    missing_text = ""
  ) %>%
  tab_style(
    style = cell_fill(color = "grey90"), 
    locations = cells_body(rows = c(1*n_models+(1:n_models), 3*n_models+(1:n_models)),  columns = everything())
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"), 
    locations = list(
      cells_body(),
      cells_column_labels(),
      cells_row_groups()
    )
  )
Sys.setenv(CHROMOTE_CHROME_ARGS = "--headless=new")
gt_table %>%
  gtsave(
    filename = "model_performance.png",
    path = "./final_outputs/figures/model_performance/",
    zoom = 3,         
    vwidth = 900,    
    vheight = 1200    
  )




# Create the table with global Serif formatting
gt_table <- gt(mp_plot) %>%
  tab_header(
    title = "Model Performance"
  ) %>%
  cols_label(
    age = "Age",
    model = "Model",
    avg_rmse = "RMSE",
    delta_rmse = "ΔRMSE",
    bmass.r2 = "R²",
    avg_auc = "AUC",
    delta_auc = "ΔAUC",
  ) %>%
  fmt_number(
    columns = c(avg_rmse, delta_rmse, bmass.r2, avg_auc, delta_auc),
    decimals = 3
  ) %>%
  tab_spanner(
    label = "Density Metrics",
    columns = c(avg_rmse, delta_rmse, bmass.r2),
  ) %>%
  tab_spanner(
    label = "Occurrence Metrics",
    columns = c(avg_auc, delta_auc),
  ) %>%
  tab_options(
    table.width = pct(60),
    column_labels.font.weight = "bold",
    data_row.padding = px(5),
    table.font.size = 16
  ) %>% 
  sub_missing(
    missing_text = ""
  ) %>%
  tab_style(
    style = cell_fill(color = "grey90"), 
    locations = cells_body(
      rows = c(1*n_models+(1:n_models), 3*n_models+(1:n_models)),  
      columns = everything()
    )
  ) %>%
  # Global Serif / Times New Roman implementation across all table elements
  opt_table_font(
    font = list(
      "Times New Roman",
      default_fonts()
    )
  )

gtsave(
  data = gt_table,
  filename = "model_performance.docx",
  path = "./final_outputs/figures/model_performance/"
)




# --------------------------------------------------------------------------------
# 7. Partial Effects Plots
# --------------------------------------------------------------------------------
# We visualize the influence of each environmental covariate on the response 
# while holding other variables constant. These plots illustrate the non-linear 
# effects of shelf distance and temperature anomalies across all hake age classes.

# Global plotting parameters for biomass shelf distance axes
bmass_shelf_dist_ylim <- c(-3, 0.5)


# 7.1 Total Hake ----------------------------------------------------------------------------


# 7.1.1 Density Model Partial Effects -------------------------------------------------------

# Extract partial effects from density model
partial_effects_bmass_total <- predict(best.gams$bmass.total, type = "terms", se.fit = TRUE)
bmass_train_dat_total <- hake_js[hake_js$wgt_total > 0,]

# --- Shelf distance partial effect ---
partial_bmass_total_sd <- data.frame(x = bmass_train_dat_total$shelf_distance_m, 
                         effect = partial_effects_bmass_total$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_bmass_total$se.fit[, "s(shelf_distance_m)"])

partial_bmass_total_sd <- partial_bmass_total_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_bmass_total_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = (x/1000)), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  scale_y_continuous(limits = bmass_shelf_dist_ylim)+
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/total/total_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/total/total_bmass_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# --- Spatially Varying (SV) Temperature Coefficient ---
# For models with spatially varying slopes, we calculate the coefficient relative to 
# the T250 anomaly and visualize areas of significant influence

slope <- as.vector(partial_effects_bmass_total$fit[,2] / bmass_train_dat_total$t250m_Ano_fill)
conf.int            <- as.vector((partial_effects_bmass_total$se.fit[,2] / bmass_train_dat_total$t250m_Ano_fill) * 1.96)
up.t                <- slope - conf.int
low.t               <- slope + conf.int
sig   <- (low.t < 0 & up.t < 0) | (low.t > 0 & up.t > 0)

vc.bmass.total <- data.frame(slope, conf.int, up.t, low.t, sig,
                             bmass_train_dat_total[c("lat","lon","easting_m", "northing_m")])

tb_t250vc <- ggplot(vc.bmass.total[vc.bmass.total$sig,]) +
  geom_sf(data = namerica_state_utm_df,
          fill = NA, color = "grey10", size = 0.2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = slope), shape = 15, size = 0.01) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  ylab("Latitude") +
  xlab("Longitude") +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = 'Longitude',
       y = 'Latitude',
       color = expression(paste('log'[e], '(B) / ºC')), 
       title = 'Total hake biomass model \n T250 variable coefficient') +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1, 2),
                        oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./final_outputs/figures/model_partial_effects/total/total_bmass_SVt250_vc.png', tb_t250vc,
       width = 2.5, height = 4.5, units = "in", bg = "white")
ggsave('./final_outputs/figures/model_partial_effects/total/total_bmass_SVt250_vc.tiff', tb_t250vc,
       device = "tiff", width = 2.5, height = 4.5, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Bivariate function of latitude and longitude spatial effect ---
slope <- as.vector(partial_effects_bmass_total$fit[,1])
vc.bmass.total <- data.frame(slope, bmass_train_dat_total[c("lat","lon","easting_m", "northing_m")])

tb_latlon <- ggplot(vc.bmass.total) +
  geom_sf(data = namerica_state_utm_df,
          fill = NA, color = "grey10", size = 0.2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = slope), shape = 15, size = 0.01) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  ylab("Latitude") +
  xlab("Longitude") +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = 'Longitude',
       y = 'Latitude',
       color = expression(paste('partial effect')), 
       title = 'Effect of s(lat,lon) on total hake density') +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1, 2),
                        oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./final_outputs/figures/model_partial_effects/total/total_bmass_latlon.png', tb_latlon,
       width = 2.5, height = 4.5, units = "in", bg = "white")
ggsave('./final_outputs/figures/model_partial_effects/total/total_bmass_latlon.tiff', tb_latlon,
       device = "tiff", width = 2.5, height = 4.5, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# 7.1.2 Occurrence Model Partial Effects -------------------------------------------------------

# Extract partial effects from Occurrence model
partial_effects_pa_total <- predict(best.gams$pa.total, type = "terms", se.fit = TRUE)

# --- Shelf Distance ---
partial_pa_total_sd <- data.frame(x = hake_js$shelf_distance_m, 
                         effect = partial_effects_pa_total$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_pa_total$se.fit[, "s(shelf_distance_m)"])
partial_pa_total_sd <- partial_pa_total_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_total_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = (x/1000)), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/total/total_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/total/total_pa_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Spatially Stationary T250 Anomaly ---
partial_pa_total_sst250 <- data.frame(x = hake_js$t250m_Ano_fill, 
                             effect = partial_effects_pa_total$fit[, "s(t250m_Ano_fill)"], 
                             se = partial_effects_pa_total$se.fit[, "s(t250m_Ano_fill)"])
partial_pa_total_sst250 <- partial_pa_total_sst250 %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_total_sst250, aes(x = x, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/total/total_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/total/total_pa_SSt250_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# 7.2 Age 2 Hake ----------------------------------------------------------------------------

# 7.2.1 Density Model Partial Effects -------------------------------------------------------

# Extract partial effects from density model
partial_effects_bmass_age2 <- predict(best.gams$bmass.age2, type = "terms", se.fit = TRUE)
bmass_train_dat_age2 <- hake_js[hake_js$age2 > 0,]

# --- Shelf Distance --- 
partial_bmass_age2_sd <- data.frame(x = bmass_train_dat_age2$shelf_distance_m, 
                         effect = partial_effects_bmass_age2$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_bmass_age2$se.fit[, "s(shelf_distance_m)"])
partial_bmass_age2_sd <- partial_bmass_age2_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_bmass_age2_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  scale_y_continuous(limits = bmass_shelf_dist_ylim)+ 
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_bmass_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Spatially Stationary T250 Anomaly ---
partial_bmass_age2_t250 <- data.frame(x = bmass_train_dat_age2$t250m_Ano_fill, 
                           effect = partial_effects_bmass_age2$fit[, "s(t250m_Ano_fill)"], 
                           se = partial_effects_bmass_age2$se.fit[, "s(t250m_Ano_fill)"])
partial_bmass_age2_t250 <- partial_bmass_age2_t250 %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_bmass_age2_t250, aes(x = x, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_bmass_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_bmass_SSt250_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# 7.2.2 Occurrence Model Partial Effects -------------------------------------------------------

# Extract partial effects from occurrence model
partial_effects_pa_age2 <- predict(best.gams$pa.age2, type = "terms", se.fit = TRUE)

# --- Shelf Distance ---
partial_pa_age2_sd <- data.frame(x = hake_js$shelf_distance_m, 
                         effect = partial_effects_pa_age2$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_pa_age2$se.fit[, "s(shelf_distance_m)"])
partial_pa_age2_sd <- partial_pa_age2_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age2_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_pa_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# --- Spatially Stationary T250 Anomaly ---
partial_pa_age2_sst250 <- data.frame(x = hake_js$t250m_Ano_fill, 
                           effect = partial_effects_pa_age2$fit[, "s(t250m_Ano_fill)"], 
                           se = partial_effects_pa_age2$se.fit[, "s(t250m_Ano_fill)"])
partial_pa_age2_sst250 <- partial_pa_age2_sst250 %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age2_sst250, aes(x = x, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age2/age2_pa_SSt250_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# 7.3 Age 3-4 Hake ----------------------------------------------------------------------------

# 7.3.1 Density Model Partial Effects -------------------------------------------------------

# Extract partial effects from density model
partial_effects_bmass_age3_4 <- predict(best.gams$bmass.age3_4, type = "terms", se.fit = TRUE)
bmass_train_dat_age3_4 <- hake_js[hake_js$age3_4 > 0,]

# --- Shelf Distance --- 
partial_bmass_age3_4_sd <- data.frame(x = bmass_train_dat_age3_4$shelf_distance_m, 
                         effect = partial_effects_bmass_age3_4$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_bmass_age3_4$se.fit[, "s(shelf_distance_m)"])

partial_bmass_age3_4_sd <- partial_bmass_age3_4_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_bmass_age3_4_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  scale_y_continuous(limits = bmass_shelf_dist_ylim)+ 
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_bmass_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# 7.3.2 Occurrence Model Partial Effects -------------------------------------------------------

# Extract partial effects from occurrence model
partial_effects_pa_age3_4 <- predict(best.gams$pa.age3_4, type = "terms", se.fit = TRUE)

# --- Shelf Distance --- 
partial_pa_age3_4_sd <- data.frame(x = hake_js$shelf_distance_m, 
                         effect = partial_effects_pa_age3_4$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_pa_age3_4$se.fit[, "s(shelf_distance_m)"])
partial_pa_age3_4_sd <- partial_pa_age3_4_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age3_4_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_pa_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Spatially Stationary T250 Anomaly ---
partial_pa_age3_4_sst250 <- data.frame(x = hake_js$t250m_Ano_fill, 
                           effect = partial_effects_pa_age3_4$fit[, "s(t250m_Ano_fill)"], 
                           se = partial_effects_pa_age3_4$se.fit[, "s(t250m_Ano_fill)"])
partial_pa_age3_4_sst250 <- partial_pa_age3_4_sst250 %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age3_4_sst250, aes(x = x, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age3_4/age3_4_pa_SSt250_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# 7.4 Age 5+ Hake ----------------------------------------------------------------------------

# 7.4.1 Density Model Partial Effects -------------------------------------------------------

# Extract partial effects from density model
partial_effects_bmass_age5_20 <- predict(best.gams$bmass.age5_20, type = "terms", se.fit = TRUE)
bmass_train_dat_age5_20 <- hake_js[hake_js$age5_20 > 0,]

# --- Shelf Distance --- 
partial_bmass_age5_20_sd <- data.frame(x = bmass_train_dat_age5_20$shelf_distance_m, 
                         effect = partial_effects_bmass_age5_20$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_bmass_age5_20$se.fit[, "s(shelf_distance_m)"])
partial_bmass_age5_20_sd <- partial_bmass_age5_20_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_bmass_age5_20_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  scale_y_continuous(limits = bmass_shelf_dist_ylim)+
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_bmass_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Spatially Varying (SV) Temperature Coefficient ---
slope <- as.vector(partial_effects_bmass_age5_20$fit[,2] / bmass_train_dat_age5_20$t250m_Ano_fill)
conf.int            <- as.vector((partial_effects_bmass_age5_20$se.fit[,2] / bmass_train_dat_age5_20$t250m_Ano_fill) * 1.96)
up.t                <- slope - conf.int
low.t               <- slope + conf.int
sig   <- (low.t < 0 & up.t < 0) | (low.t > 0 & up.t > 0)

vc.bmass.age5_20 <- data.frame(slope, conf.int, up.t, low.t, sig,
                             bmass_train_dat_age5_20[c("lat","lon","easting_m", "northing_m")])

a5_20_t250vc <- ggplot(vc.bmass.age5_20[vc.bmass.age5_20$sig, ]) +
  geom_sf(data = namerica_state_utm_df,
          fill = NA, color = "grey10", size = 0.2) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = slope), shape = 15, size = 0.01) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
  ylab("Latitude") +
  xlab("Longitude") +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  labs(x = 'Longitude',
       y = 'Latitude',
       color = expression(paste('log'[e], '(B) / ºC')), 
       title = 'Age 5-20 hake biomass model \n T250 variable coefficient') +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1, 2),
                        oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
        panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./final_outputs/figures/model_partial_effects/age5_20/age5_20_bmass_SVt250_vc.png', a5_20_t250vc,
       width = 2.5, height = 4.5, units = "in", bg = "white")
ggsave('./final_outputs/figures/model_partial_effects/age5_20/age5_20_bmass_SVt250_vc.tiff', a5_20_t250vc,
       device = "tiff", width = 2.5, height = 4.5, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# 7.4.2 Occurrence Model Partial Effects -------------------------------------------------------

# Extract partial effects from Occurrence model
partial_effects_pa_age5_20 <- predict(best.gams$pa.age5_20, type = "terms", se.fit = TRUE)

# --- Shelf Distance --- 
partial_pa_age5_20_sd <- data.frame(x = hake_js$shelf_distance_m, 
                         effect = partial_effects_pa_age5_20$fit[, "s(shelf_distance_m)"], 
                         se = partial_effects_pa_age5_20$se.fit[, "s(shelf_distance_m)"])
partial_pa_age5_20_sd <- partial_pa_age5_20_sd %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age5_20_sd, aes(x = x/1000, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x/1000), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of Shelf Distance",
       x = "km from 200 meter shelf-break",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_pa_shelf_dist_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Spatially Stationary T250 Anomaly ---
partial_pa_age5_20_sst250 <- data.frame(x = hake_js$t250m_Ano_fill, 
                           effect = partial_effects_pa_age5_20$fit[, "s(t250m_Ano_fill)"], 
                           se = partial_effects_pa_age5_20$se.fit[, "s(t250m_Ano_fill)"])
partial_pa_age5_20_sst250 <- partial_pa_age5_20_sst250 %>%
  mutate(lower = effect - 1.96 * se,
         upper = effect + 1.96 * se)

p <- ggplot(partial_pa_age5_20_sst250, aes(x = x, y = effect)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+
  geom_line(color = "black") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "grey", alpha = 0.7) +
  geom_rug(aes(x = x), sides = "b", color = "black", linewidth = 0.1) + # Add rug plot
  labs(title = "Partial Effect of t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(text = element_text(family = "serif"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")
ggsave("./final_outputs/figures/model_partial_effects/age5_20/age5_20_pa_SSt250_PE.tiff", p,
       device = "tiff", width = 4, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")



# --------------------------------------------------------------------------------
# 8. Spatially Explicit Forecast Maps
# --------------------------------------------------------------------------------
# We visualize the mean biomass forecasts, prediction uncertainty (CV), 
# and projected biomass changes relative to the historical period


# 8.1: 2024 Mean Biomass Forecast ------------------------------------------------
# Shows a map of the predicted biomass (kg) for August 2024 
# across our study region for each age group
o1 <-  ggplot(spatExpl_sim_out) +
  geom_sf(data = namerica_state_utm_df,
          fill = "grey90", color = "grey10", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = mean_post), shape = 15, size = 0.06) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  geom_line(data = border, aes(x = easting_m, y = northing_m), color = 'red')+
  facet_wrap(~ age, nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       color = 'kg',
       title = paste0('Mean of hake biomass prediction simulations for August 2024')) +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  scale_color_viridis(limits = mean_lims, oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
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

ggsave('./final_outputs/figures/forecasts/mean_pred_2024.png', o1,
       width = 4.5, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/mean_pred_2024.tiff', o1,
       device = "tiff", width = 4.5, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

  
# 8.2: Forecast Uncertainty (CV) -----------------------------------------------
# The Coefficient of Variation (CV) identifies the precision of 
# our Monte Carlo simulations across the study area with respect to the prediction magnitude (CV = SE/Mean)
o2 <- ggplot(spatExpl_sim_out) +
         geom_sf(data = namerica_state_utm_df,
            fill = "grey90", color = "grey10", size = 0.2) +
         geom_point(aes(x = easting_m, y = northing_m, color = cv_post), shape = 15, size = 0.06) +
         geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
         facet_wrap(~ age, nrow = 1) +
         labs(x = "Longitude",
              y = "Latitude",
              color = 'SE',
              title = paste0('CV of hake biomass prediction simulations for August 2024')) +
         scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
         ylim(ymin, ymax) +
         scale_color_viridis(limits = cv_lims, oob = scales::squish) +
         theme_simple(base_size = 9, grid = TRUE) +
         theme(text = element_text(family = "serif"),
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

ggsave('./final_outputs/figures/forecasts/CV_2024.png', o2,
       width = 5.5, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/CV_2024.tiff', o2,
       device = "tiff", width = 5.5, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# 8.3: Forecast Uncertainty (SE) -----------------------------------------------
# The Standard Error (SE) identifies the precision of 
# our Monte Carlo simulations across the study area 
o3 <- ggplot(spatExpl_sim_out) +
  geom_sf(data = namerica_state_utm_df,
          fill = "grey90", color = "grey10", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = se_post), shape = 15, size = 0.06) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  facet_wrap(~ age, nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       color = 'SE',
       title = paste0('SE of hake biomass prediction simulations for August 2024')) +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  scale_color_gradient2(limits = c(0,3e4), oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
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

ggsave('./final_outputs/figures/forecasts/SE_2024.png', o3,
       width = 5.5, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/SE_2024.tiff', o3,
       device = "tiff", width = 5.5, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# 8.4: Projected Distribution Shifts -------------------------------------------
# Visualizing the net change in biomass (kg) between pre_year and post_year
o4 <- ggplot(spatExpl_sim_out) +
  geom_sf(data = namerica_state_utm_df,
          fill = "grey90", color = "grey10", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = mean_diff), shape = 15, size = 0.06) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  facet_wrap(~ age, nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       color = 'kg',
       title = paste0('Difference between predictions of 2023 and 2024 biomass')) +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-4e3,4e3),
                        midpoint = 0,
                        oob = scales::squish) +    
  theme_simple(base_size = 9, grid = TRUE) +
  theme(text = element_text(family = "serif"),
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
ggsave('./final_outputs/figures/forecasts/diff_hist_fut.png', o4,
       width = 5.5, height = 3, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/diff_hist_fut.tiff', o4,
       device = "tiff", width = 5.5, height = 3, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

# --------------------------------------------------------------------------------
# 9. Cohort-Specific Forecast Profiles 
# --------------------------------------------------------------------------------
# We generate multi-panel figures for each age class to showcase the predicted 
# biomass in August 2024, model uncertainty (CV) of that prediction, and net biomass change in 2024
# relative to the 2023 baseline

# --- Total Hake Profile ---
total_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                          mean_lims = c(0,120),
                                          se_lims = c(0,60), 
                                          change_lims = c(-40, 40),
                                          # change_lims = c(-60, 60), 
                                          age_group = "Total")
ggsave('./final_outputs/figures/forecasts/age_class_maps/total_ageclass_maps.png', total_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/age_class_maps/total_ageclass_maps.tiff', total_ageclass_maps,
       device = "tiff", width = 5.7, height = 2.8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

# --- Age 2 Hake Profile ---
age2_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                         mean_lims = c(0,6),
                                         se_lims = c(0,3), 
                                         change_lims = c(-2, 2),
                                         age_group = "Age 2")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age2_ageclass_maps.png', age2_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age2_ageclass_maps.tiff', age2_ageclass_maps,
       device = "tiff", width = 5.7, height = 2.8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

# --- Age 3-4 Hake Profile ---
age3_4_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                           mean_lims = c(0,50),
                                           se_lims = c(0,25), 
                                           change_lims = c(-15, 15),
                                           age_group = "Age 3-4")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age3_4_ageclass_maps.png', age3_4_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age3_4_ageclass_maps.tiff', age3_4_ageclass_maps,
       device = "tiff", width = 5.7, height = 2.8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")


# --- Age 5+ Hake Profile ---
age5_20_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                            mean_lims = c(0,60),
                                            se_lims = c(0,30), 
                                            change_lims = c(-20,20),
                                            age_group = "Age 5+")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age5_20_ageclass_maps.png', age5_20_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/age_class_maps/age5_20_ageclass_maps.tiff', age5_20_ageclass_maps,
       device = "tiff", width = 5.7, height = 2.8, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")




plot_year_forecast <- function(plot_year, forecast_output = spat_expl_simulation_outputs, lims){
  
  plot_2015 <- ggplot(forecast_output %>% filter(forecast_output$age == "Total" & year == plot_year)) +
    geom_point(aes(x = easting_m, y = northing_m, color = mean/1000), shape = 15, size = 0.06) +
    labs(x = "Longitude",
         y = "Latitude",
         color = 'kg x 10^3',
         title = paste0(plot_year)) +
    geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +  
    geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", size = 0.2) +
    scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
    ylim(ymin, ymax) +
    scale_color_viridis(limits = lims, oob = scales::squish) +
    theme_simple(base_size = 8, grid = TRUE) +
    theme(text = element_text(family = "serif"),
          plot.title = element_text(hjust = 0.5),
          plot.title.position = "panel",
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
  
}


p_2015 <- plot_year_forecast(plot_year = 2015, lims = c(0,120))
p_2017 <- plot_year_forecast(plot_year = 2017, lims = c(0,120))
p_2019 <- plot_year_forecast(plot_year = 2019, lims = c(0,120))


ggsave('./final_outputs/figures/forecasts/2015_total_forecast.tiff', p_2015,
       device = "tiff", width = 1.4, height = 3.75, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

ggsave('./final_outputs/figures/forecasts/2017_total_forecast.tiff', p_2017,
       device = "tiff", width = 1.4, height = 3.75, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

ggsave('./final_outputs/figures/forecasts/2019_total_forecast.tiff', p_2019,
       device = "tiff", width = 1.4, height = 3.75, units = "in", 
       dpi = 500, compression = "lzw", bg = "white")

  









# 
# # showcase difference between 2017 and 2019
# total_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
#                                           mean_lims = c(0,120),
#                                           se_lims = c(0,60), 
#                                           change_lims = c(-60, 60), 
#                                           age_group = "Total")
# age2_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
#                                          mean_lims = c(0,1),
#                                          se_lims = c(0,0.5), 
#                                          change_lims = c(-0.5,0.5),
#                                          age_group = "Age 2")
# age3_4_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
#                                            mean_lims = c(0,20),
#                                            se_lims = c(0,10), 
#                                            change_lims = c(-10, 10),
#                                            age_group = "Age 3-4")
# age5_20_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
#                                             mean_lims = c(0,60),
#                                             se_lims = c(0,30), 
#                                             change_lims = c(-30,30),
#                                             age_group = "Age 5+")
# ggsave('C:/Users/jakem/Downloads/down.png', age5_20_ageclass_maps,
#        width = 5.7, height = 2.8, units = "in", bg = "white")
# ggsave('C:/Users/jakem/Downloads/down.png', age5_20_ageclass_maps,
#        device = "tiff", width = 5.7, height = 2.8, units = "in", 
#        dpi = 500, compression = "lzw", bg = "white")







# --------------------------------------------------------------------------------
# 10. Regional Aggregation & Transboundary Analysis 
# --------------------------------------------------------------------------------
# We summarize the Delta Model simulation results into aggregate biomass signals 
# to compare predicted shifts between the US and Canadian EEZs.


# --- Aggregated Biomass in U.S. and Canada --- 
coeff = 1/10
o5 <- ggplot(regional_sim_out, aes(x = year)) +
  geom_vline(xintercept = c(2015, 2019), alpha = 0.5) +
  annotate('rect', xmin = 2015, xmax = 2019, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.3) + 
  geom_line(data = regional_sim_out[regional_sim_out$region == "US",], aes(y = mean_reg, color = "US"), linetype = "solid", linewidth = 0.6) +
  geom_line(data = regional_sim_out[regional_sim_out$region == "CAN",], aes(y = mean_reg/coeff, color = "CAN"), linetype = "solid", linewidth = 0.6) +
  geom_line(data = hake_obs_reg[hake_obs_reg$region == "US",], aes(x = year, y = wgts, color = "US"), linetype = "dashed", linewidth = 0.6) +
  geom_line(data = hake_obs_reg[hake_obs_reg$region == "CAN",], aes(x = year, y = wgts, color = "CAN"), linetype = "dashed", linewidth = 0.6) +
  facet_wrap(~ age, nrow = 2) +
  scale_x_continuous(breaks = regional_sim_out$year) +  # Set breaks for x-axis
  scale_y_continuous(
    name = "kg in the U.S.",
    sec.axis = sec_axis(~ .*coeff, name = "kg in Canada")
  ) +
  labs(title = 'Predicted aggregated hake biomass signal in the US and Canadian EEZs',
       x = 'Year Observed') +
  scale_color_manual(name = "Region", values = c("US" = "blue", "CAN" = "red")) +
  theme_minimal() +
  theme(text = element_text(family = "serif"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y.right = element_text(margin = margin(l = 5), color = "red"), 
        axis.text.y.left = element_text(margin = margin(l = 5), color = "blue"),
        panel.spacing.x = unit(1, "lines"))  
ggsave('./final_outputs/figures/forecasts/regional_biomass.png', o5,
       width = 7, height = 5, units = "in", bg = "white")
ggsave('./final_outputs/figures/forecasts/regional_biomass.tiff', o5,
       device = "tiff", width = 7, height = 5, units = "in",
       dpi = 500, compression = "lzw", bg = "white")


# end of script

