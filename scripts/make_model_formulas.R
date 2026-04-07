#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Model Formula Definition & Generation
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Defines a suite of candidate Generalized Additive Model (GAM) 
#'              formulas for model selection. This script generates formulas 
#'              for both occurrence (presence/absence) and density (log-biomass) 
#'              across four age classes: Total, Age 2, Age 3-4, and Age 5-20 (also called Age 5+).
#' Output: ./gam/formulas.RData
#' #############################################################################

if(!dir.exists("./gam/")) dir.create("./gam/")

# ----------------------------------------------------------------------------
# 1. Make model formulas for Total biomass 
# ----------------------------------------------------------------------------
formula_strings <- c( 
  
  'log_wgt_total ~ s(easting_m, northing_m,  k = 25) +
                   s(shelf_distance_m, k = 4) +
                   s(year_fac, bs = "re")',
  
  'log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                   s(t100m_Ano_fill, k = 5) +
                   s(shelf_distance_m, k = 4) +
                   s(year_fac, bs = "re")',
  
  'log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                   s(t250m_Ano_fill, k = 5) +
                   s(shelf_distance_m, k = 4) +
                   s(year_fac, bs = "re")',
  
  'log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                   s(easting_m, northing_m, by = t100m_Ano_fill, k = 10) +
                   s(shelf_distance_m, k = 4) +
                   s(year_fac, bs = "re")',
  
  'log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                   s(easting_m, northing_m, by = t250m_Ano_fill, k = 10) +
                   s(shelf_distance_m, k = 4) +
                   s(year_fac, bs = "re")' 
)

formulas.bmass.total <- list()

for(i in 1:length(formula_strings)){
  formulas.bmass.total[[i]] <- formula(formula_strings[i])
}

formulas.pa.total <- lapply(formulas.bmass.total, FUN = function(x) update(x, present_total ~ .))

# ----------------------------------------------------------------------------
# 2. Age 2
# ----------------------------------------------------------------------------
formulas.bmass.age2 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age2 ~ .))
formulas.pa.age2 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age2 ~ .))

# ----------------------------------------------------------------------------
# 3. Age 3-4
# ----------------------------------------------------------------------------
formulas.bmass.age3_4 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age3_4 ~ .))
formulas.pa.age3_4 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age3_4 ~ .))

# ----------------------------------------------------------------------------
# 4. Age 5
# ----------------------------------------------------------------------------
formulas.bmass.age5_20 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age5_20 ~ .))
formulas.pa.age5_20 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age5_20 ~ .))

# ----------------------------------------------------------------------------
# 5. Save formulas
# ----------------------------------------------------------------------------
formulas <- list(
  bmass.total = formulas.bmass.total,
  pa.total = formulas.pa.total,
  bmass.age2 = formulas.bmass.age2, 
  pa.age2 = formulas.pa.age2,
  bmass.age3_4 = formulas.bmass.age3_4, 
  pa.age3_4 = formulas.pa.age3_4,
  bmass.age5_20 = formulas.bmass.age5_20,
  pa.age5_20 = formulas.pa.age5_20
)

save(formulas, file = "./gam/formulas.RData")

# end of script


