#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Model Performance Evaluation & Selection (LFO Cross-Validation)
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Performs Leave-Future-Out (LFO) cross-validation over a 7-year
#'              retrospective period to evaluate biomass (RMSE) and 
#'              occurrence (AUC) models. Identifies the best-performing 
#'              GAMs for each age class and saves the fully parameterized 
#'              versions for final forecasting.
#' Output: ./outputs/model_performance/performance_7year_array.RData
#'         ./gam/best.gams.RData
#' #############################################################################

# ----------------------------------------------------------------------------
# 1. Configuration and Data Initialization
# ----------------------------------------------------------------------------

num_years <- 7

# Define output subdirectories
dirs <- c("./outputs/model_performance/total", 
          "./outputs/model_performance/age2", 
          "./outputs/model_performance/age3_4", 
          "./outputs/model_performance/age5_20")
invisible(lapply(dirs, function(d) if(!dir.exists(d)) dir.create(d, recursive = TRUE)))

# Load required assets
load("./gam/formulas.RData")


# ----------------------------------------------------------------------------
# 2. Retrospective Density Evaluation
# ----------------------------------------------------------------------------
# Evaluate density models by training on n years and calculating RMSE 
# against observed biomass in year n+1.


# --- Total Biomass ---
RMSE_js_total_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                                      target = 'log_wgt_total', 
                                                      subset = 'wgt_total', 
                                                      fxn = train_js_bmass_gams, 
                                                      forms = formulas$bmass.total,
                                                      nyears = num_years)

save(RMSE_js_total_array_7years, file = './outputs/model_performance/total/RMSE_js_total_array_7years.RData')

# --- Age 2 ---
RMSE_js_age2_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                                     target = 'log_age2', 
                                                     subset = 'age2', 
                                                     fxn = train_js_bmass_gams, 
                                                     forms = formulas$bmass.age2, 
                                                     nyears = num_years)

save(RMSE_js_age2_array_7years, file = './outputs/model_performance/age2/RMSE_js_age2_array_7years.RData')

# --- Age 3-4 ---
RMSE_js_age3_4_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                                       target = 'log_age3_4', 
                                                       subset = 'age3_4', 
                                                       fxn = train_js_bmass_gams, 
                                                       forms = formulas$bmass.age3_4,
                                                       nyears = num_years)

save(RMSE_js_age3_4_array_7years, file = './outputs/model_performance/age3_4/RMSE_js_age3_4_array_7years.RData')

# --- Age 5+ ---
RMSE_js_age5_20_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                                        target = 'log_age5_20', 
                                                        subset = 'age5_20', 
                                                        fxn = train_js_bmass_gams, 
                                                        forms = formulas$bmass.age5_20,
                                                        nyears = num_years)

save(RMSE_js_age5_20_array_7years, file = './outputs/model_performance/age5_20/RMSE_js_age5_20_array_7years.RData')


# ----------------------------------------------------------------------------
# 3. Retrospective Occurrence Evaluation
# ----------------------------------------------------------------------------
# Evaluate occurrence models by training on n years and calculating AUC 
# against observed occurrence in year n+1.

# --- Total --- 
AUC_js_total_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                                  target = 'present_total', 
                                                  fxn = train_js_pa_gams, 
                                                  forms = formulas$pa.total,
                                                  nyears = num_years)

save(AUC_js_total_array_7years, file = './outputs/model_performance/total/AUC_js_total_array_7years.RData')

# --- Age 2 --- 
AUC_js_age2_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                                 target = 'present_age2', 
                                                 fxn = train_js_pa_gams, 
                                                 forms = formulas$pa.age2,
                                                 nyears = num_years)

save(AUC_js_age2_array_7years, file = './outputs/model_performance/age2/AUC_js_age2_array_7years.RData')

# --- Age 3-4 --- 
AUC_js_age3_4_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                                   target = 'present_age3_4', 
                                                   fxn = train_js_pa_gams, 
                                                   forms = formulas$pa.age3_4,
                                                   nyears = num_years)

save(AUC_js_age3_4_array_7years, file = './outputs/model_performance/age3_4/AUC_js_age3_4_array_7years.RData')

# --- Age 5+ --- 
AUC_js_age5_20_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                                    target = 'present_age5_20', 
                                                    fxn = train_js_pa_gams, 
                                                    forms = formulas$pa.age5_20,
                                                    nyears = num_years)

save(AUC_js_age5_20_array_7years, file = './outputs/model_performance/age5_20/AUC_js_age5_20_array_7years.RData')


# Organize outputs
performance_7year_array <- list()
performance_7year_array$total_bmass <- RMSE_js_total_array_7years
performance_7year_array$age2_bmass <- RMSE_js_age2_array_7years
performance_7year_array$age3_4_bmass <- RMSE_js_age3_4_array_7years
performance_7year_array$age5_20_bmass <- RMSE_js_age5_20_array_7years
performance_7year_array$total_pa <- AUC_js_total_array_7years
performance_7year_array$age2_pa <- AUC_js_age2_array_7years
performance_7year_array$age3_4_pa <- AUC_js_age3_4_array_7years
performance_7year_array$age5_20_pa <- AUC_js_age5_20_array_7years
save(performance_7year_array, file = './outputs/model_performance/performance_7year_array.RData')

# ----------------------------------------------------------------------------
# 4. Model Selection
# ----------------------------------------------------------------------------
# Determine best models based on lowest RMSE/highest AUC in above LFO analyses

# Find best models based on lowest RMSE over the retrospective period
ind_best <- list()
ind_best$total_bmass <- which(RMSE_js_total_array_7years[,ncol(RMSE_js_total_array_7years)] == min(RMSE_js_total_array_7years[,ncol(RMSE_js_total_array_7years)]))
ind_best$age2_bmass <- which(RMSE_js_age2_array_7years[,ncol(RMSE_js_age2_array_7years)] == min(RMSE_js_age2_array_7years[,ncol(RMSE_js_total_array_7years)]))
ind_best$age3_4_bmass <- which(RMSE_js_age3_4_array_7years[,ncol(RMSE_js_age3_4_array_7years)] == min(RMSE_js_age3_4_array_7years[,ncol(RMSE_js_age3_4_array_7years)]))
ind_best$age5_20_bmass <- which(RMSE_js_age5_20_array_7years[,ncol(RMSE_js_age5_20_array_7years)] == min(RMSE_js_age5_20_array_7years[,ncol(RMSE_js_age5_20_array_7years)]))

# Find best pa models based on highest AUC over the retrospective period
ind_best$total_pa <- which(AUC_js_total_array_7years[,ncol(AUC_js_total_array_7years)] == max(AUC_js_total_array_7years[,ncol(AUC_js_total_array_7years)]))
ind_best$age2_pa <- which(AUC_js_age2_array_7years[,ncol(AUC_js_age2_array_7years)] == max(AUC_js_age2_array_7years[,ncol(AUC_js_age2_array_7years)]))
ind_best$age3_4_pa <- which(AUC_js_age3_4_array_7years[,ncol(AUC_js_age3_4_array_7years)] == max(AUC_js_age3_4_array_7years[,ncol(AUC_js_age3_4_array_7years)]))
ind_best$age5_20_pa <- which(AUC_js_age5_20_array_7years[,ncol(AUC_js_age5_20_array_7years)] == max(AUC_js_age5_20_array_7years[,ncol(AUC_js_age5_20_array_7years)]))


# ----------------------------------------------------------------------------
# 5. Save Formulas
# ----------------------------------------------------------------------------

best.formulas <- list()
best.formulas[['bmass.total']] <- formulas$bmass.total[[ind_best$total_bmass]]
best.formulas[['bmass.age2']] <- formulas$bmass.age2[[ind_best$age2_bmass]]
best.formulas[['bmass.age3_4']] <- formulas$bmass.age3_4[[ind_best$age3_4_bmass]] ## problem -- best model doesn't take temp into consideration
best.formulas[['bmass.age5_20']] <- formulas$bmass.age5_20[[ind_best$age5_20_bmass]]

best.formulas[['pa.total']] <- formulas$pa.total[[ind_best$total_pa]]
best.formulas[['pa.age2']] <- formulas$pa.age2[[ind_best$age2_pa]]
best.formulas[['pa.age3_4']] <- formulas$pa.age3_4[[ind_best$age3_4_pa]]
best.formulas[['pa.age5_20']] <- formulas$pa.age5_20[[ind_best$age5_20_pa]]

save(best.formulas, file = './gam/best.formulas.RData')


# ----------------------------------------------------------------------------
# 6. Final Model Training
# ----------------------------------------------------------------------------
# Parameterizing best models with all available data for final forecasting.

models <-  paste0("Model ", 1:length(formulas$bmass.total))

# --- Total Density --- 
fully_trained_total_bmass <- train_js_bmass_gams(data = hake_js, subset = 'wgt_total', formulas = formulas$bmass.total)
total_bmass_perform <- data.frame(aic = double(length = length(fully_trained_total_bmass)), 
                                  r2 = double(length = length(fully_trained_total_bmass)),
                                  dev.expl= double(length = length(fully_trained_total_bmass)),
                                  model = models,
                                  age = "Total", 
                                  avg_rmse = performance_7year_array$total_bmass[,ncol(performance_7year_array$total_bmass)])

for(i in 1:length(fully_trained_total_bmass)){
  total_bmass_perform$aic[i] <- AIC(fully_trained_total_bmass[[i]])
  total_bmass_perform$r2[i] <- summary(fully_trained_total_bmass[[i]])$r.sq
  total_bmass_perform$dev.expl[i] <- summary(fully_trained_total_bmass[[i]])$dev.expl
}

# --- Total Occurrence --- 
fully_trained_total_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.total)
total_pa_perform <- data.frame(aic = double(length = length(fully_trained_total_pa)), 
                               model = models, 
                               age = "Total",
                               avg_auc = performance_7year_array$total_pa[,ncol(performance_7year_array$total_pa)])
for(i in 1:length(fully_trained_total_pa)){
  total_pa_perform$aic[i] <- AIC(fully_trained_total_bmass[[i]])
}


# --- Age 2 Density --- 
fully_trained_age2_bmass <- train_js_bmass_gams(data = hake_js, subset = 'age2', formulas = formulas$bmass.age2)
age2_bmass_perform <- data.frame(aic = double(length = length(fully_trained_age2_bmass)), 
                                 r2 = double(length = length(fully_trained_age2_bmass)),
                                 dev.expl= double(length = length(fully_trained_age2_bmass)),
                                 model = models,
                                 age = "Age 2",
                                 avg_rmse = performance_7year_array$age2_bmass[,ncol(performance_7year_array$age2_bmass)])
for(i in 1:length(fully_trained_age2_bmass)){
  age2_bmass_perform$aic[i] <- AIC(fully_trained_age2_bmass[[i]])
  age2_bmass_perform$r2[i] <- summary(fully_trained_age2_bmass[[i]])$r.sq
  age2_bmass_perform$dev.expl[i] <- summary(fully_trained_age2_bmass[[i]])$dev.expl
}

# --- Age 2 Occurrence --- 
fully_trained_age2_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age2)
age2_pa_perform <- data.frame(aic = double(length = length(fully_trained_age2_pa)), 
                              model = models,
                              age = "Age 2",
                              avg_auc = performance_7year_array$age2_pa[,ncol(performance_7year_array$age2_pa)])
for(i in 1:length(fully_trained_age2_pa)){
  age2_pa_perform$aic[i] <- AIC(fully_trained_age2_pa[[i]])
}


# --- Age 3 Density --- 
fully_trained_age3_4_bmass <- train_js_bmass_gams(data = hake_js, subset = 'age3_4', formulas = formulas$bmass.age3_4)
age3_4_bmass_perform <- data.frame(aic = double(length = length(fully_trained_age3_4_bmass)), 
                                   r2 = double(length = length(fully_trained_age3_4_bmass)),
                                   dev.expl= double(length = length(fully_trained_age3_4_bmass)),
                                   model = models,
                                   age = "Age 3-4",
                                   avg_rmse = performance_7year_array$age3_4_bmass[,ncol(performance_7year_array$age3_4_bmass)])
for(i in 1:length(fully_trained_age3_4_bmass)){
  age3_4_bmass_perform$aic[i] <- AIC(fully_trained_age3_4_bmass[[i]])
  age3_4_bmass_perform$r2[i] <- summary(fully_trained_age3_4_bmass[[i]])$r.sq
  age3_4_bmass_perform$dev.expl[i] <- summary(fully_trained_age3_4_bmass[[i]])$dev.expl
}

# --- Age 3 Occurrence --- 
fully_trained_age3_4_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age3_4)
age3_4_pa_perform <- data.frame(aic = double(length = length(fully_trained_age3_4_pa)), 
                                model = models,
                                age = "Age 3-4",
                                avg_auc = performance_7year_array$age3_4_pa[,ncol(performance_7year_array$age3_4_pa)])
for(i in 1:length(fully_trained_age3_4_pa)){
  age3_4_pa_perform$aic[i] <- AIC(fully_trained_age3_4_pa[[i]])
}


# --- Age 5+ Density --- 
fully_trained_age5_20_bmass <- train_js_bmass_gams(data = hake_js, subset = 'age5_20', formulas = formulas$bmass.age5_20)
age5_20_bmass_perform <- data.frame(aic = double(length = length(fully_trained_age5_20_bmass)), 
                                    r2 = double(length = length(fully_trained_age5_20_bmass)),
                                    dev.expl= double(length = length(fully_trained_age5_20_bmass)),
                                    model = models,
                                    age = "Age 5+",
                                    avg_rmse = performance_7year_array$age5_20_bmass[,ncol(performance_7year_array$age5_20_bmass)])
for(i in 1:length(fully_trained_age5_20_bmass)){
  age5_20_bmass_perform$aic[i] <- AIC(fully_trained_age5_20_bmass[[i]])
  age5_20_bmass_perform$r2[i] <- summary(fully_trained_age5_20_bmass[[i]])$r.sq
  age5_20_bmass_perform$dev.expl[i] <- summary(fully_trained_age5_20_bmass[[i]])$dev.expl
}

# --- Age 5+ Occurrence --- 
fully_trained_age5_20_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age5_20)
age5_20_pa_perform <- data.frame(aic = double(length = length(fully_trained_age5_20_pa)), 
                                 model = models,
                                 age = "Age 5+",
                                 avg_auc = performance_7year_array$age5_20_pa[,ncol(performance_7year_array$age5_20_pa)])
for(i in 1:length(fully_trained_age5_20_pa)){
  age5_20_pa_perform$aic[i] <- AIC(fully_trained_age5_20_pa[[i]])
}


# ----------------------------------------------------------------------------
# 7. Save
# ----------------------------------------------------------------------------

# Save performance metrics of best models
model.performance <- list()
model.performance$total.bmass <- total_bmass_perform
model.performance$age2.bmass <- age2_bmass_perform
model.performance$age3_4.bmass <- age3_4_bmass_perform
model.performance$age5_20.bmass <- age5_20_bmass_perform

model.performance$total.pa <- total_pa_perform
model.performance$age2.pa <- age2_pa_perform
model.performance$age3_4.pa <- age3_4_pa_perform
model.performance$age5_20.pa <- age5_20_pa_perform

save(model.performance, file = './outputs/model_performance/model.performance.RData')


# Save best models
best.gams <- list()
best.gams[['bmass.total']] <- fully_trained_total_bmass[[ind_best$total_bmass]]
best.gams[['bmass.age2']] <- fully_trained_age2_bmass[[ind_best$age2_bmass]]
best.gams[['bmass.age3_4']] <- fully_trained_age3_4_bmass[[ind_best$age3_4_bmass]] ## problem -- best model doesn't take temp into consideration
best.gams[['bmass.age5_20']] <- fully_trained_age5_20_bmass[[ind_best$age5_20_bmass]]

best.gams[['pa.total']] <- fully_trained_total_pa[[ind_best$total_pa]]
best.gams[['pa.age2']] <- fully_trained_age2_pa[[ind_best$age2_pa]]
best.gams[['pa.age3_4']] <- fully_trained_age3_4_pa[[ind_best$age3_4_pa]]
best.gams[['pa.age5_20']] <- fully_trained_age5_20_pa[[ind_best$age5_20_pa]]

save(best.gams, file = './gam/best.gams.RData')


# end of script
