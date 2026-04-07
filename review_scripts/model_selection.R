
##############################################################################
# This script analyzes model performance of biomass and presence models separately
# using a leave-future-out cross-validation procedure.
# Best models for each age-class are designated based on performance during the
# past n years of predictions 
# (set as n = 7 to analyze models across warm and cold years in past)
num_years <- 7


# make directories
if(!dir.exists("./outputs/"))
{dir.create("./outputs/")}
if(!dir.exists("./outputs/model_performance/"))
{dir.create("./outputs/model_performance/")}


# load model formulas to be tested
load("./gam/formulas.RData")


############################################################################
############################################################################
# Biomass Model Performance
############################################################################
############################################################################

RMSE_js_total_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                              target = 'log_wgt_total', 
                                              subset = 'wgt_total', 
                                              fxn = train_js_bmass_gams, 
                                              forms = formulas$bmass.total,
                                              nyears = num_years)

if(!dir.exists("./outputs/model_performance/total/"))
  {dir.create("./outputs/model_performance/total/")}

save(RMSE_js_total_array_7years, file = './outputs/model_performance/total/RMSE_js_total_array_7years.RData')


RMSE_js_age2_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                             target = 'log_age2', 
                                             subset = 'age2', 
                                             fxn = train_js_bmass_gams, 
                                             forms = formulas$bmass.age2, 
                                             nyears = num_years)

if(!dir.exists("./outputs/model_performance/age2/"))
  {dir.create("./outputs/model_performance/age2/")}

save(RMSE_js_age2_array_7years, file = './outputs/model_performance/age2/RMSE_js_age2_array_7years.RData')


RMSE_js_age3_4_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                               target = 'log_age3_4', 
                                               subset = 'age3_4', 
                                               fxn = train_js_bmass_gams, 
                                               forms = formulas$bmass.age3_4,
                                               nyears = num_years)

if(!dir.exists("./outputs/model_performance/age3_4/"))
  {dir.create("./outputs/model_performance/age3_4/")}

save(RMSE_js_age3_4_array_7years, file = './outputs/model_performance/age3_4/RMSE_js_age3_4_array_7years.RData')


RMSE_js_age5_20_array_7years <- lfo_analysis_bmass_gams(data = hake_js, 
                                                target = 'log_age5_20', 
                                                subset = 'age5_20', 
                                                fxn = train_js_bmass_gams, 
                                                forms = formulas$bmass.age5_20,
                                                nyears = num_years)

if(!dir.exists("./outputs/model_performance/age5_20/"))
  {dir.create("./outputs/model_performance/age5_20/")}

save(RMSE_js_age5_20_array_7years, file = './outputs/model_performance/age5_20/RMSE_js_age5_20_array_7years.RData')



############################################################################
############################################################################
# Presence Model Performance
############################################################################
############################################################################

AUC_js_total_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                            target = 'present_total', 
                                            fxn = train_js_pa_gams, 
                                            forms = formulas$pa.total,
                                            nyears = num_years)

save(AUC_js_total_array_7years, file = './outputs/model_performance/total/AUC_js_total_array_7years.RData')


AUC_js_age2_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                           target = 'present_age2', 
                                           fxn = train_js_pa_gams, 
                                           forms = formulas$pa.age2,
                                           nyears = num_years)

save(AUC_js_age2_array_7years, file = './outputs/model_performance/age2/AUC_js_age2_array_7years.RData')


AUC_js_age3_4_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                             target = 'present_age3_4', 
                                             fxn = train_js_pa_gams, 
                                             forms = formulas$pa.age3_4,
                                             nyears = num_years)

save(AUC_js_age3_4_array_7years, file = './outputs/model_performance/age3_4/AUC_js_age3_4_array_7years.RData')


AUC_js_age5_20_array_7years <- lfo_analysis_pa_gams(data = hake_js, 
                                            target = 'present_age5_20', 
                                            fxn = train_js_pa_gams, 
                                            forms = formulas$pa.age5_20,
                                            nyears = num_years)

save(AUC_js_age5_20_array_7years, file = './outputs/model_performance/age5_20/AUC_js_age5_20_array_7years.RData')


# load above 
# load('./outputs/model_performance/age5_20/AUC_js_age5_20_array_7years.RData')
# load('./outputs/model_performance/age5_20/RMSE_js_age5_20_array_7years.RData')
# load('./outputs/model_performance/age3_4/AUC_js_age3_4_array_7years.RData')
# load('./outputs/model_performance/age3_4/RMSE_js_age3_4_array_7years.RData')
# load('./outputs/model_performance/age2/AUC_js_age2_array_7years.RData')
# load('./outputs/model_performance/age2/RMSE_js_age2_array_7years.RData')
# load('./outputs/model_performance/total/AUC_js_total_array_7years.RData')
# load('./outputs/model_performance/total/RMSE_js_total_array_7years.RData')

# Organize into a more manageable form
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



########################################################################
########################################################################
## Determine best models based on lowest RMSE/highest AUC in above
## LFO analyses
########################################################################

# find best models based on lowest RMSE over the retrospective period
ind_best <- list()
ind_best$total_bmass <- which(RMSE_js_total_array_7years[,ncol(RMSE_js_total_array_7years)] == min(RMSE_js_total_array_7years[,ncol(RMSE_js_total_array_7years)]))
ind_best$age2_bmass <- which(RMSE_js_age2_array_7years[,ncol(RMSE_js_age2_array_7years)] == min(RMSE_js_age2_array_7years[,ncol(RMSE_js_total_array_7years)]))
ind_best$age3_4_bmass <- which(RMSE_js_age3_4_array_7years[,ncol(RMSE_js_age3_4_array_7years)] == min(RMSE_js_age3_4_array_7years[,ncol(RMSE_js_age3_4_array_7years)]))
ind_best$age5_20_bmass <- which(RMSE_js_age5_20_array_7years[,ncol(RMSE_js_age5_20_array_7years)] == min(RMSE_js_age5_20_array_7years[,ncol(RMSE_js_age5_20_array_7years)]))

# find best pa models based on highest AUC over the retrospective period
ind_best$total_pa <- which(AUC_js_total_array_7years[,ncol(AUC_js_total_array_7years)] == max(AUC_js_total_array_7years[,ncol(AUC_js_total_array_7years)]))
ind_best$age2_pa <- which(AUC_js_age2_array_7years[,ncol(AUC_js_age2_array_7years)] == max(AUC_js_age2_array_7years[,ncol(AUC_js_age2_array_7years)]))
ind_best$age3_4_pa <- which(AUC_js_age3_4_array_7years[,ncol(AUC_js_age3_4_array_7years)] == max(AUC_js_age3_4_array_7years[,ncol(AUC_js_age3_4_array_7years)]))
ind_best$age5_20_pa <- which(AUC_js_age5_20_array_7years[,ncol(AUC_js_age5_20_array_7years)] == max(AUC_js_age5_20_array_7years[,ncol(AUC_js_age5_20_array_7years)]))



########################################################################
## Save best formulas
########################################################################
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



########################################################################
########################################################################
## Parameterize best models with all training data and find performance metrics (R2, dev explained, AIC)
## Save best models for later use
#################################################################

models <-  paste0("Model ", 1:length(formulas$bmass.total))

# total hake
#################################################################
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

fully_trained_total_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.total)
total_pa_perform <- data.frame(aic = double(length = length(fully_trained_total_pa)), 
                               model = models, 
                               age = "Total",
                               avg_auc = performance_7year_array$total_pa[,ncol(performance_7year_array$total_pa)])
for(i in 1:length(fully_trained_total_pa)){
  total_pa_perform$aic[i] <- AIC(fully_trained_total_bmass[[i]])
}


# age 2
#################################################################
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

fully_trained_age2_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age2)
age2_pa_perform <- data.frame(aic = double(length = length(fully_trained_age2_pa)), 
                              model = models,
                              age = "Age 2",
                              avg_auc = performance_7year_array$age2_pa[,ncol(performance_7year_array$age2_pa)])
for(i in 1:length(fully_trained_age2_pa)){
  age2_pa_perform$aic[i] <- AIC(fully_trained_age2_pa[[i]])
}


# age 3
######################################################
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

fully_trained_age3_4_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age3_4)
age3_4_pa_perform <- data.frame(aic = double(length = length(fully_trained_age3_4_pa)), 
                                model = models,
                                age = "Age 3-4",
                                avg_auc = performance_7year_array$age3_4_pa[,ncol(performance_7year_array$age3_4_pa)])
for(i in 1:length(fully_trained_age3_4_pa)){
  age3_4_pa_perform$aic[i] <- AIC(fully_trained_age3_4_pa[[i]])
}


# age 5-20
######################################################
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

fully_trained_age5_20_pa <- train_js_pa_gams(data = hake_js, formulas = formulas$pa.age5_20)
age5_20_pa_perform <- data.frame(aic = double(length = length(fully_trained_age5_20_pa)), 
                                 model = models,
                                 age = "Age 5+",
                                 avg_auc = performance_7year_array$age5_20_pa[,ncol(performance_7year_array$age5_20_pa)])
for(i in 1:length(fully_trained_age5_20_pa)){
  age5_20_pa_perform$aic[i] <- AIC(fully_trained_age5_20_pa[[i]])
}


# save performance metrics of best models
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


# save best models
########################################################################
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
