
## Make figures for publishing

## Uses output from 
## - match_data.R
## - model_selection.R
## - generate_output.R
####################################################################

## Sections
###########
## #   Load files and set important plotting variables 
## # 0 Location Map 
## # 1 Hake Data 
## # 2 Temperature Forecasts 
## # 3 Model performance (RMSE and AUC)
## # 4 Partial Effects Plots
## # 5 Output maps
## # 6 Regionally aggregated output


####################################################################
####################################################################
## Load files and set important plotting variables 
####################################################################
####################################################################
setwd("C:/Users/jakem/OneDrive/Desktop/GitHub/hake_j-scope_sdms/")

if(!dir.exists('./outputs/figures/'))
  {dir.create('./outputs/figures/')}

# load files and make vars for plotting
load('./data/hake_js.RData')
load('./outputs/output_stats/spat_expl_simulation_outputs.RData')
load('./outputs/output_stats/spatExpl_sim_out.RData')
load('./outputs/output_stats/region_diff.RData')
load('./outputs/output_stats/regional_sim_out.RData')
load('./gam/best.gams.RData')


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



####################################################################
####################################################################
## 0 Show geographical location with respect to US West Coast
## (match_data.R)
####################################################################
####################################################################

o1 <- ggplot() +
  geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
          linewidth = 0.2) +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax= ymax), 
            color = "red", fill = NA) +
  # scale_x_continuous(breaks = c(seq(-200, -50, 15)), limits = c(xmin-1000000, xmax+1000000)) +
  # scale_y_continuous(breaks = c(seq(20, 70, 10)), limits = c(ymin-1500000, ymax+1000000)) +
  scale_x_continuous(breaks = c(seq(-200, -50, 15)), limits = c(xmin-500000*2, xmax+700000*2)) +
  scale_y_continuous(breaks = c(seq(20, 70, 10)), limits = c(ymin-1000000*2, ymax+1000000*1.5)) +
  labs(x = "Longitude",
       y = "Latitude") +
  # theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/placeofstudy.png', o1,
       width = 3, height = 4, units = "in", bg = "white")



####################################################################
####################################################################
## # 1 Hake Data 
## (match_data.R)
####################################################################
####################################################################

if(!dir.exists("./outputs/figures/data/"))
  {dir.create("./outputs/figures/data/")}


# plot weight over all years for each age class
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
       color = expression(paste('log','(B)'))) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

# # old color scheme -- trash
# h1 <- ggplot(hake_dat) +
#   geom_point(aes(x = easting_m, y = northing_m, color = log_bmass), shape = 15, size = 0.2) +
#   geom_sf(data = namerica_state_utm_df, fill = "grey90", color = "grey30",
#           size = 0.2) +
#   geom_sf(data = bath_200_utm_df, fill = NA, color = "grey20", size = 0.2) +
#   scale_color_viridis(na.value = 'transparent', limits = c(0,17), oob = scales::squish) +
#   facet_wrap( ~ age, nrow = 1) +
#   scale_x_continuous(breaks = c(-128, -126, -124), limits = c(xmin, xmax)) +
#   ylim(ymin, ymax) +
#   labs(x = "Longitude",
#        y = "Latitude",
#        color = expression(paste('log'[e], '(B)'))) +
#   theme_simple(base_size = 9, grid = TRUE) +
#   theme(panel.spacing.x = unit(-0.5, "pt"),
#         panel.spacing = unit(1, "lines"),
#         legend.justification = c(0, 0),
#         legend.position = c(0.01, 0.005),
#         legend.title.align = 0.5,
#         legend.key.size = unit(0.4, "cm"),
#         legend.background = element_blank(),
#         legend.text = element_text(size = 6),
#         legend.title = element_text(size = 6),
#         panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/data/hake_by_ageclass.png', h1,
       width = 7, height = 3, units = "in", bg = "white")


# hake by year
# plotted latitude against distance to shelf
h1_3 <- ggplot(hake_js) +
  geom_hline(aes(yintercept = 48.25), linetype = "dashed", alpha = 0.7) +
  geom_vline(aes(xintercept = 0), linetype = "solid", alpha = 0.7) +
  geom_point(aes(x = shelf_distance_m/1000, y = lat, color = log_wgt_total), size = 0.00001) +
  scale_color_viridis(na.value = "transparent") +
  scale_x_reverse() +
  facet_wrap(~year_fac, nrow = 4) +
  labs(x = "Shelf Distance (km)",
       y = "Latitude (ºN)", 
       color = expression(paste('log'[e], '(B)')))+
  theme_simple(base_size = 12, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 7),
        legend.title = element_text(size = 7),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/data/hake_lat_shelfdist.png', h1_3,
       width = 6, height = 8, units = "in", bg = "white")


####################################################################
####################################################################
## # 2 Temperature Forecasts 
## (match_data.R)
####################################################################
####################################################################

# ready the data for plotting
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

# plot t250 over all years
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
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/data/js_t250_all_years.png', t1,
       width = 7, height = 4.5, units = "in", bg = "white")

# plot t100 over all years
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
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/data/js_t100_all_years.png', t2,
       width = 7, height = 4.5, units = "in", bg = "white")

####################################################################
####################################################################
# plot hot year, cold year, and 2024 temperature
####################################################################
####################################################################

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
cold_ind <- order(mean_temps$t250, decreasing = F)[1] # coldest year spatially averaged
hot_ind <- order(mean_temps$t250, decreasing = F)[nrow(mean_temps)] # hottest year spatially averaged
fut_ind <- which(mean_temps$year == 2024) 
years_hcf <- mean_temps$year[c(hot_ind, cold_ind, fut_ind)]

dat_hcf <- hake_dat[is.element(hake_dat$year_fac, c(years_hcf)), ]
dat_hcf$year_fac <- factor(dat_hcf$year_fac, levels = c(years_hcf)) # data.frame with only hot, neutral, and cold years


# plot t250 over hot/cold/2024 years
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
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/data/js_t250_hot_cold_fut_years.png', t4,
       width = 4, height = 3, units = "in", bg = "white")

# plot t100 over hot/cold/2024 years
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
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/data/js_t100_hot_cold_fut_years.png', t5,
       width = 4, height = 3, units = "in", bg = "white")


####################################################################
####################################################################
## # 3 LFO-CV Model Performance Heatmaps (RMSE and AUC)
## (model_selection.R)
####################################################################
####################################################################
if(!dir.exists('./outputs/figures/model_performance/'))
  {dir.create('./outputs/figures/model_performance/')}
if(!dir.exists('./outputs/figures/model_performance/total/'))
  {dir.create('./outputs/figures/model_performance/total/')}
if(!dir.exists('./outputs/figures/model_performance/age2/'))
  {dir.create('./outputs/figures/model_performance/age2/')}
if(!dir.exists('./outputs/figures/model_performance/age3_4/'))
  {dir.create('./outputs/figures/model_performance/age3_4/')}
if(!dir.exists('./outputs/figures/model_performance/age5_20/'))
  {dir.create('./outputs/figures/model_performance/age5_20/')}

load('./outputs/model_performance/performance_7year_array.RData')

#make plotting variables and load formulas
n_models <- nrow(performance_7year_array$total_bmass)
model_names <- paste0("Model ", n_models:1) # c("Model 5", "Model 4", "Model 3", "Model 2", "Model 1") 
pred_years_7years <- c(tail(unique(hake_js$year), 7), 'mean')

rmse_lims <- c(1,3)
auc_lims <- c(0.5,1)

# total 
####################################################################
# biomass

# last 7 years biomass heatmap
mp_t_1 <- plot_heatmap(data = performance_7year_array$total_bmass, years = pred_years_7years, 
                       names = model_names, type = 'RMSE', lims = rmse_lims, 
                       age = "Total")
ggsave("./outputs/figures/model_performance/total/total_bmass_RMSE_7years.png", mp_t_1,
       width = 6, height = 3, units = "in", bg = "white")

# last 7 years presence heatmap
mp_t_2 <- plot_heatmap(data = performance_7year_array$total_pa, years = pred_years_7years, 
                       names = model_names, type = 'AUC', lims = auc_lims, 
                       age = "Total")
ggsave("./outputs/figures/model_performance/total/total_AUC_7years.png", mp_t_2,
       width = 6, height = 3, units = "in", bg = "white")


# age 2
####################################################################
# biomass

# last 7 years biomass heatmap
mp_a2_1 <- plot_heatmap(data = performance_7year_array$age2_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 2")
ggsave("./outputs/figures/model_performance/age2/age2_bmass_RMSE_7years.png", mp_a2_1,
       width = 6, height = 3, units = "in", bg = "white")

# last 7 years presence heatmap
mp_a2_2 <- plot_heatmap(data = performance_7year_array$age2_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 2")
ggsave("./outputs/figures/model_performance/age2/age2_AUC_7years.png", mp_a2_2,
       width = 6, height = 3, units = "in", bg = "white")


# age 3-4
####################################################################
# biomass

# last 7 years biomass heatmap
mp_a3_1 <- plot_heatmap(data = performance_7year_array$age3_4_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 3-4")
ggsave("./outputs/figures/model_performance/age3_4/age3_4_bmass_RMSE_7years.png", mp_a3_1,
       width = 6, height = 3, units = "in", bg = "white")

# last 7 years presence heatmap
mp_a3_2 <- plot_heatmap(data = performance_7year_array$age3_4_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 3-4")
ggsave("./outputs/figures/model_performance/age3_4/age3_4_AUC_7years.png", mp_a3_2,
       width = 6, height = 3, units = "in", bg = "white")


# age 5-20
####################################################################
# biomass

# last 7 years biomass heatmap
mp_a5_1 <- plot_heatmap(data = performance_7year_array$age5_20_bmass, years = pred_years_7years, 
                        names = model_names, type = 'RMSE', lims = rmse_lims, 
                        age = "Age 5+")
ggsave("./outputs/figures/model_performance/age5_20/age5_20_bmass_RMSE_7years.png", mp_a5_1,
       width = 6, height = 3, units = "in", bg = "white")

# last 7 years presence heatmap
mp_a5_2 <- plot_heatmap(data = performance_7year_array$age5_20_pa, years = pred_years_7years, 
                        names = model_names, type = 'AUC', lims = auc_lims, 
                        age = "Age 5+")
ggsave("./outputs/figures/model_performance/age5_20/age5_20_AUC_7years.png", mp_a5_2,
       width = 6, height = 3, units = "in", bg = "white")



# fully trained model performance
####################################################################
load('./outputs/model_performance/model.performance.RData')



for(i in 1:length(model.performance)){
  model.performance[[i]] <- model.performance[[i]][1:9,]
}


## Reorder and condense metrics for plotting
## put presence metrics in same model rows as biomass metrics
mp <- bind_rows(model.performance)

# mp$pa.aic <- c(mp$aic[which(!is.na(mp$avg_auc))], rep(NA, 24))
mp$avg_auc <- mp$avg_auc[!is.na(mp$avg_auc)]
names(mp) <- c("bmass.aic", "bmass.r2", "bmass.dev.expl", "model", "age", "avg_rmse", "avg_auc") #, "pa.aic")
mp$age <- ifelse((1:nrow(mp)) %% n_models == 1, mp$age, NA)
mp <- mp %>% filter(!is.na(avg_rmse))

## add delta RMSE/AUC
################################################################
# delta rmse
delta_rmse <- double(length = nrow(mp))

delta_rmse[1:n_models] <- mp$avg_rmse[1:n_models] - min(mp$avg_rmse[1:n_models])
delta_rmse[1*n_models+(1:n_models)] <- mp$avg_rmse[1*n_models+(1:n_models)] - min(mp$avg_rmse[1*n_models+(1:n_models)])
delta_rmse[2*n_models+(1:n_models)] <- mp$avg_rmse[2*n_models+(1:n_models)] - min(mp$avg_rmse[2*n_models+(1:n_models)])
delta_rmse[3*n_models+(1:n_models)] <- mp$avg_rmse[3*n_models+(1:n_models)] - min(mp$avg_rmse[3*n_models+(1:n_models)])

# delta auc
delta_auc <- double(length = nrow(mp))
delta_auc[1:n_models] <- mp$avg_auc[1:n_models] - max(mp$avg_auc[1:n_models])  
delta_auc[1*n_models+(1:n_models)] <- mp$avg_auc[1*n_models+(1:n_models)] - max(mp$avg_auc[1*n_models+(1:n_models)]) 
delta_auc[2*n_models+(1:n_models)] <-  mp$avg_auc[2*n_models+(1:n_models)] - max(mp$avg_auc[2*n_models+(1:n_models)])
delta_auc[3*n_models+(1:n_models)] <- mp$avg_auc[3*n_models+(1:n_models)] - max(mp$avg_auc[3*n_models+(1:n_models)])

mp$delta_rmse <- delta_rmse
mp$delta_auc <- delta_auc
mp_plot <- mp[, c("age", "model", "avg_rmse", "delta_rmse", "bmass.r2", "avg_auc","delta_auc")]

# # only for presentation
# #######################################
# mp_plot$model[which(mp_plot$model == 'Model 1')] <- 'Spatial Only'
# mp_plot$model[which(mp_plot$model == 'Model 2')] <- 'S.S. T100'
# mp_plot$model[which(mp_plot$model == 'Model 3')] <- 'S.S. T250'
# mp_plot$model[which(mp_plot$model == 'Model 4')] <- 'S.V. T100'
# mp_plot$model[which(mp_plot$model == 'Model 5')] <- 'S.V. T250'
# ###############################

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
  # tab_style(
  #   style = list(cell_text(weight = "bold")),#, cell_fill(color = "lightcyan")),
  #   locations = cells_body(rows = c(5, 8, 11, 20), columns = c("model","avg_rmse", "delta_rmse", "bmass.r2"))
  # ) %>%
  # tab_style(
  #   style = list(cell_text(weight = "bold")),#, cell_fill(color = "#F9E3D6")),
  #   locations = cells_body(rows = c(3, 8, 13, 18), columns = c("model","avg_auc", "delta_auc"))
  # ) %>%
  # tab_style(
  #   style = list(cell_text(weight = "bold")),#, cell_fill(color = "violet")),
  #   locations = cells_body(rows = c(8), columns = c("model"))
  # ) %>% 
  tab_style(
    style = cell_fill(color = "grey90"), 
    locations = cells_body(rows = c(1*n_models+(1:n_models), 3*n_models+(1:n_models)),  columns = everything())
  ) %>%
  tab_style(
    style = cell_text(font = "Times New Roman"),  # Specify the font family
    locations = list(
      cells_body(),
      cells_column_labels(),
      cells_row_groups()
    )
  )
gtsave(gt_table, file = './outputs/figures/model_performance/model_performance.html')
webshot(url = './outputs/figures/model_performance/model_performance.html',
        file = './outputs/figures/model_performance/model_performance.png',
        vwidth = 800, vheight = 800)

####################################################################
####################################################################
## # 4 Partial Effects Plots
# - using only selected models from model_selection.R
# (model_selection.R)
####################################################################
####################################################################

load('./gam/best.gams.RData')


if(!dir.exists("./outputs/figures/model_partial_effects/"))
{dir.create("./outputs/figures/model_partial_effects/")}
if(!dir.exists("./outputs/figures/model_partial_effects/total/"))
{dir.create("./outputs/figures/model_partial_effects/total/")}
if(!dir.exists("./outputs/figures/model_partial_effects/age2/"))
{dir.create("./outputs/figures/model_partial_effects/age2/")}
if(!dir.exists("./outputs/figures/model_partial_effects/age3_4/"))
{dir.create("./outputs/figures/model_partial_effects/age3_4/")}
if(!dir.exists("./outputs/figures/model_partial_effects/age5_20/"))
{dir.create("./outputs/figures/model_partial_effects/age5_20/")}

bmass_shelf_dist_ylim <- c(-3, 0.5)

##################################################
## Total Hake
##################################################

# biomass
##################################################

# Extract partial effects
partial_effects_bmass_total <- predict(best.gams$bmass.total, type = "terms", se.fit = TRUE)
bmass_train_dat_total <- hake_js[hake_js$wgt_total > 0,]

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/total/total_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# S.V. Temperature biomass model
###################################################
# Get variable coefficient for t250 in biomass model
slope <- as.vector(partial_effects_bmass_total$fit[,2] / bmass_train_dat_total$t250m_Ano_fill)
conf.int            <- as.vector((partial_effects_bmass_total$se.fit[,2] / bmass_train_dat_total$t250m_Ano_fill) * 1.96)
up.t                <- slope - conf.int
low.t               <- slope + conf.int
sig   <- (low.t < 0 & up.t < 0) | (low.t > 0 & up.t > 0)

vc.bmass.total <- data.frame(slope, conf.int, up.t, low.t, sig,
                             bmass_train_dat_total[c("lat","lon","easting_m", "northing_m")])


## Visualize temperature variable coefficient
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/model_partial_effects/total/total_bmass_SVt250_vc.png', tb_t250vc,
       width = 2.5, height = 4.5, units = "in", bg = "white")


# visualize s(lat, lon) spatial effect
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/model_partial_effects/total/total_bmass_latlon.png', tb_latlon,
       width = 2.5, height = 4.5, units = "in", bg = "white")



# presence 
##################################################

# Extract partial effects
partial_effects_pa_total <- predict(best.gams$pa.total, type = "terms", se.fit = TRUE)

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/total/total_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")

# SS t250
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
  labs(title = "Partial Effect of SS t250",
       x = "ºC",
       y = "Partial Effect") +
  theme_classic(base_size = 10)+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/total/total_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


##################################################
## age 2
##################################################

# biomass
##################################################

# Extract partial effects
partial_effects_bmass_age2 <- predict(best.gams$bmass.age2, type = "terms", se.fit = TRUE)
bmass_train_dat_age2 <- hake_js[hake_js$age2 > 0,]

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age2/age2_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")

# SS t250
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age2/age2_bmass_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# presence 
##################################################

# Extract partial effects
partial_effects_pa_age2 <- predict(best.gams$pa.age2, type = "terms", se.fit = TRUE)

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age2/age2_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# SS t250
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age2/age2_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


##################################################
## age 3-4
##################################################

# biomass
##################################################

# Extract partial effects
partial_effects_bmass_age3_4 <- predict(best.gams$bmass.age3_4, type = "terms", se.fit = TRUE)
bmass_train_dat_age3_4 <- hake_js[hake_js$age3_4 > 0,]

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age3_4/age3_4_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")

# presence 
##################################################

# Extract partial effects
partial_effects_pa_age3_4 <- predict(best.gams$pa.age3_4, type = "terms", se.fit = TRUE)

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age3_4/age3_4_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# SS t250
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age3_4/age3_4_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


##################################################
## age 5-20
##################################################

# biomass
##################################################

# Extract partial effects
partial_effects_bmass_age5_20 <- predict(best.gams$bmass.age5_20, type = "terms", se.fit = TRUE)
bmass_train_dat_age5_20 <- hake_js[hake_js$age5_20 > 0,]

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age5_20/age5_20_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


## Get variable coefficient for t250 in biomass model
slope <- as.vector(partial_effects_bmass_age5_20$fit[,2] / bmass_train_dat_age5_20$t250m_Ano_fill)
conf.int            <- as.vector((partial_effects_bmass_age5_20$se.fit[,2] / bmass_train_dat_age5_20$t250m_Ano_fill) * 1.96)
up.t                <- slope - conf.int
low.t               <- slope + conf.int
sig   <- (low.t < 0 & up.t < 0) | (low.t > 0 & up.t > 0)

vc.bmass.age5_20 <- data.frame(slope, conf.int, up.t, low.t, sig,
                             bmass_train_dat_age5_20[c("lat","lon","easting_m", "northing_m")])

#Visualize results
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/model_partial_effects/age5_20/age5_20_bmass_SVt250_vc.png', a5_20_t250vc,
       width = 2.5, height = 4.5, units = "in", bg = "white")


# presence 
##################################################

# Extract partial effects
partial_effects_pa_age5_20 <- predict(best.gams$pa.age5_20, type = "terms", se.fit = TRUE)

# shelf dist
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age5_20/age5_20_pa_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# SS t250
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
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/age5_20/age5_20_pa_SSt250_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")





####################################################################
####################################################################
## # 5 Output Maps (Spatially Explicit Output)
# - mean simulated delta prediction for 2024
# - SE and standardized SE of simulated delta prediction for 2024
# - straight delta predictions
# (generate_output.R)
####################################################################
####################################################################

if(!dir.exists("./outputs/figures/forecasts/"))
  {dir.create("./outputs/figures/forecasts/")}

load('./outputs/output_stats/spatExpl_sim_out.RData')
load('./outputs/output_stats/hist_fut_temp.RData')



cv_lims <- c(0, 0.5)
mean_lims <- c(0, 3e4)
delta_lims <- c(0, 5e4)

# mean simulated delta prediction for 2024 across age class
o1 <-  ggplot(spatExpl_sim_out) +
  geom_sf(data = namerica_state_utm_df,
          fill = "grey90", color = "grey10", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = sim_delta_2024), shape = 15, size = 0.06) +
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/forecasts/mean_pred_2024.png', o1,
       width = 4.5, height = 3, units = "in", bg = "white")

  
# CV of simulated delta prediction for 2024
o2 <- ggplot(spatExpl_sim_out) +
         geom_sf(data = namerica_state_utm_df,
            fill = "grey90", color = "grey10", size = 0.2) +
         geom_point(aes(x = easting_m, y = northing_m, color = sim_cv_2024), shape = 15, size = 0.06) +
         geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
         facet_wrap(~ age, nrow = 1) +
         labs(x = "Longitude",
              y = "Latitude",
              color = 'SE',
              title = paste0('SE of hake biomass prediction simulations for August 2024')) +
         scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
         ylim(ymin, ymax) +
         scale_color_viridis(limits = cv_lims, oob = scales::squish) +
         theme_simple(base_size = 9, grid = TRUE) +
         theme(panel.spacing.x = unit(-0.5, "pt"),
              panel.spacing = unit(1, "lines"),
              legend.justification = c(0, 0),
              legend.position = c(0.01, 0.005),
              legend.title.align = 0.5,
              legend.key.size = unit(0.4, "cm"),
              legend.background = element_blank(),
              legend.text = element_text(size = 6),
              legend.title = element_text(size = 6),
              panel.grid.major = element_line(size = 0.2))

ggsave('./outputs/figures/forecasts/CV_2024.png', o2,
       width = 5.5, height = 3, units = "in", bg = "white")

# change from 2023 predictions to 2024
o4 <- ggplot(spatExpl_sim_out) +
  geom_sf(data = namerica_state_utm_df,
          fill = "grey90", color = "grey10", size = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = delta_change), shape = 15, size = 0.06) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "grey90", size = 0.2) +
  facet_wrap(~ age, nrow = 1) +
  labs(x = "Longitude",
       y = "Latitude",
       color = 'kg',
       title = paste0('Difference between 2015 - 2023 predicted and 2024 predicted biomass')) +
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/forecasts/diff_hist_fut.png', o4,
       width = 5.5, height = 3, units = "in", bg = "white")




## Plot forecasted biomass for warm years (2015, 2017, 2019)
spat_expl_simulation_outputs$easting_m <- spatExpl_sim_out$easting_m
spat_expl_simulation_outputs$northing_m <- spatExpl_sim_out$northing_m

hist_total <- plot_hist_years_delta(data = spat_expl_simulation_outputs, 
                                    age = 'total', 
                                    delta_lims = c(0,12e4), 
                                    years = c(2015, 2017, 2019))

hist_age2 <- plot_hist_years_delta(data = spat_expl_simulation_outputs, 
                                    age = 'age 2', 
                                    delta_lims = c(0,6e3), 
                                    years = c(2015, 2017, 2019))


hist_age3_4 <- plot_hist_years_delta(data = spat_expl_simulation_outputs, 
                                   age = 'age 3-4', 
                                   delta_lims = c(0,5e4), 
                                   years = c(2015, 2017, 2019))


hist_age5_20 <- plot_hist_years_delta(data = spat_expl_simulation_outputs, 
                                   age = 'age 5-20', 
                                   delta_lims = c(0,6e4), 
                                   years = c(2015, 2017, 2019))


gridExtra::grid.arrange(grobs = hist_total, nrow = 1)
gridExtra::grid.arrange(grobs = hist_age2,nrow = 1)
gridExtra::grid.arrange(grobs = hist_age3_4,nrow = 1)
gridExtra::grid.arrange(grobs = hist_age5_20,nrow = 1)




## Plot forecasted biomass, CV, and biomass change from 2023 to 2024 for each age class
###################################################################

if(!dir.exists("./outputs/figures/forecasts/age_class_maps/"))
  {dir.create("./outputs/figures/forecasts/age_class_maps/")}

total_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                           delta_lims = c(0,120),
                                           cv_lims = cv_lims, 
                                           change_lims = c(-60, 60), 
                                           age = "total")
ggsave('./outputs/figures/forecasts/age_class_maps/total_ageclass_maps.png', total_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")

age2_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                          delta_lims = c(0,6),
                                          cv_lims = cv_lims, 
                                          change_lims = c(-3, 3), 
                                          age = "age 2")
ggsave('./outputs/figures/forecasts/age_class_maps/age2_ageclass_maps.png', age2_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")

age3_4_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                            delta_lims = c(0,50),
                                            cv_lims = cv_lims, 
                                            change_lims = c(-25, 25), 
                                          age = "age 3-4")
ggsave('./outputs/figures/forecasts/age_class_maps/age3_4_ageclass_maps.png', age3_4_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")

age5_20_ageclass_maps <- plot_ageclass_maps(data = spatExpl_sim_out, 
                                             delta_lims = c(0,60),
                                             cv_lims = cv_lims, 
                                             change_lims = c(-30,30), 
                                             age = "age 5-20")
ggsave('./outputs/figures/forecasts/age_class_maps/age5_20_ageclass_maps.png', age5_20_ageclass_maps,
       width = 5.7, height = 2.8, units = "in", bg = "white")


## find difference between 2023 and 2024 t250
hist_fut_temp_diff <- ggplot(hist_fut_temp) +
  geom_sf(data = namerica_state_utm_df,fill = "grey90", color = "grey10", linewidth = 0.2) +
  geom_point(aes(x = easting_m, y = northing_m, color = t250_change), shape = 15, size = 0.06) +
  geom_sf(data = bath_200_utm_df, fill = NA, color = "black", size = 0.2) +  
  labs(x = "Longitude",
       y = "Latitude",
       color = 'ºC',
       title = paste0('Predicted Biomass    Prediction CV   Change in Biomass')) +
  scale_x_continuous(breaks = c(-127, -125), limits = c(xmin, xmax)) +
  ylim(ymin, ymax) +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-0.5, 0.5),
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
ggsave('./outputs/figures/forecasts/hist_fut_temp_diff.png', hist_fut_temp_diff,
       width = 2, height = 2.8, units = "in", bg = "white")




####################################################################
## # 6 Regionally aggregated signal
# - sum mean simulated delta prediction using LFYO over US and CAN EEZs
# (generate_output.R)
####################################################################
load('./outputs/output_stats/regional_sim_out.RData')
load('./outputs/output_stats/region_diff.RData')
load('./outputs/output_stats/hake_obs_reg_diff.RData')
load('./outputs/output_stats/hake_obs_reg.RData')

# both observations and predictions of U.S.-Canada ratio
o7_ratio <- ggplot() +
  geom_vline(xintercept = c(2015, 2019), alpha = 0.5) +
  annotate('rect', xmin = 2015, xmax = 2019, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) + 
  geom_line(data = region_diff, aes(x = year, y = ratio_US, color = "Predicted"), linetype = "solid", linewidth = 1) +
  geom_line(data = hake_obs_reg_diff, aes(x = year, y = ratio_US, color = "Observed"), linetype = "solid", linewidth = 1) +
  facet_wrap(~ age, nrow = 2) +
  scale_x_continuous(breaks = unique(region_diff$year)) +  # Set breaks for x-axis
  scale_y_continuous(name = "hake ratio",) +
  labs(title = 'Predicted and observed ratio hake in US',
       x = 'Year Observed') +  
  scale_color_manual(name = "", values = c("Predicted" = "deepskyblue3", "Observed" = "black")) +#c("Predicted" = "deepskyblue3", "Observed" = "azure4")) +
  theme_minimal(base_size = 10)+
  theme(strip.text = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y.right = element_text(margin = margin(l = 5)),
        panel.spacing.x = unit(1.5, "lines")) # Increase left margin of right axis text
ggsave('./outputs/figures/forecasts/pred_v_obs_ratio.png', o7_ratio,
       width = 7, height = 4, units = "in", bg = "white")

# ratio of hake observations made in the U.S. to the total observations in the study area
o7_US_obs_ratio <- ggplot() +
  geom_vline(xintercept = c(2015, 2019), alpha = 0.5) +
  annotate('rect', xmin = 2015, xmax = 2019, ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.2) + 
  geom_line(data = hake_obs_reg, aes(x = year, y = ratio_num_obs_US), linetype = "solid", linewidth = 1) +
  scale_x_continuous(breaks = unique(hake_obs_reg$year)) +  # Set breaks for x-axis
  scale_y_continuous(limits = c(0.125, 1),name = "observation ratio",) +
  labs(title = 'Ratio of hake survey observations in the US',
       x = 'Year Observed') +  
  theme_minimal(base_size = 10)+
  theme(strip.text = element_text(size = 11),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y.right = element_text(margin = margin(l = 5))) # Increase left margin of right axis text
ggsave('./outputs/figures/forecasts/US_obs_ratio.png', o7_US_obs_ratio,
       width = 5, height = 3, units = "in", bg = "white")

## plot US and CAN predictions on same plot
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
    name = "kg in US",
    sec.axis = sec_axis(~ .*coeff, name = "kg in Canada")
  ) +
  labs(title = 'Predicted aggregated hake biomass signal in the US and Canadian EEZs',
       x = 'Year Observed') +
  scale_color_manual(name = "Region", values = c("US" = "blue", "CAN" = "red")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y.right = element_text(margin = margin(l = 5), color = "red"),  # Color right y-axis labels red
        axis.text.y.left = element_text(margin = margin(l = 5), color = "blue"))  # Color right y-axis labels red
ggsave('./outputs/figures/forecasts/regional_biomass.png', o5,
       width = 7, height = 5, units = "in", bg = "white")




# end of script























#######
# looking into alternative gam effects for NOAA Review


##################################################
## Total Hake
##################################################

# Extract partial effects
partial_effects_bmass_total <- predict(best.gams$bmass.total, type = "terms", se.fit = TRUE)
bmass_train_dat_total <- hake_js[hake_js$wgt_total > 0,]

# shelf dist
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
  # scale_y_continuous(limits = c(ymin, ymax))+
  theme_classic(base_size = 10)+
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(size = 12))
ggsave("./outputs/figures/model_partial_effects/total/total_bmass_shelf_dist_PE.png", p,
       width = 4, height = 3, units = "in", bg = "white")


# shelf-latitude-varying Temperature effect in biomass model
###################################################

partial_effects_bmass_total <- predict(best.gams$bmass.total, type = "terms", se.fit = TRUE)
bmass_train_dat_total <- hake_js[hake_js$wgt_total > 0,]

# Get variable coefficient for t250 in biomass model
slope <- as.vector(partial_effects_bmass_total$fit[,"s(shelf_distance_m,northing_m):t250m_Ano_fill"] / bmass_train_dat_total$t250m_Ano_fill)
conf.int            <- as.vector((partial_effects_bmass_total$se.fit[,"s(shelf_distance_m,northing_m):t250m_Ano_fill"] / bmass_train_dat_total$t250m_Ano_fill) * 1.96)
up.t                <- slope - conf.int
low.t               <- slope + conf.int
sig   <- (low.t < 0 & up.t < 0) | (low.t > 0 & up.t > 0)

vc.bmass.total <- data.frame(slope, conf.int, up.t, low.t, sig,
                             bmass_train_dat_total[c("lat","lon","easting_m", "northing_m")])


## Visualize temperature variable coefficient
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
       title = 'Total hake biomass model \n shelf/lat varying T250 coefficient') +
  scale_color_gradient2(low = "steelblue",
                        mid = "grey90",
                        high = "red3",
                        space = "Lab",
                        na.value = "transparent",
                        limits = c(-1, 2),
                        oob = scales::squish) +
  theme_simple(base_size = 9, grid = TRUE) +
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
# ggsave('./outputs/figures/model_partial_effects/total/total_bmass_SVt250_vc.png', tb_t250vc,
#        width = 2.5, height = 4.5, units = "in", bg = "white")


# visualize s(lat, lon) spatial effect
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
  theme(panel.spacing.x = unit(-0.5, "pt"),
        panel.spacing = unit(1, "lines"),
        legend.justification = c(0, 0),
        legend.position = c(0.01, 0.005),
        legend.title.align = 0.5,
        legend.key.size = unit(0.4, "cm"),
        legend.background = element_blank(),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 9),
        panel.grid.major = element_line(size = 0.2))
ggsave('./outputs/figures/model_partial_effects/total/total_bmass_latlon.png', tb_latlon,
       width = 2.5, height = 4.5, units = "in", bg = "white")
















