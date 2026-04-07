############################################################################
## Make model formulas to be used in model selection
## - each of 5 different models for presence and biomass for each age class group
############################################################################

# Make model formulas for Total biomass -- to be used for all ages models for density and occurrence
formula_strings <- c( 
                      
                      # Spatial Only
                      "log_wgt_total ~ s(easting_m, northing_m,  k = 25) +
                                       s(shelf_distance_m, k = 4) +
                                       s(year_fac, bs = 're')",
                      
                      # New -- Spatial Only, only shelf dist and northing
                      # just latxshelfdist, shelf dist
                      "log_wgt_total ~ s(shelf_distance_m, northing_m, k = 25) +
                                       s(year_fac, bs = 're')",
                      
                      # # latxlon, S.S. t100, shelf dist
                      # "log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                      #                  s(t100m_Ano_fill, k = 5) +
                      #                  s(shelf_distance_m, k = 4) +
                      #                  s(year_fac, bs = 're')",

                      
                      # S.S. T250
                      "log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                                       s(t250m_Ano_fill, k = 5) +
                                       s(shelf_distance_m, k = 4) +
                                       s(year_fac, bs = 're')",
                      
                      # New -- S.S. T250 w/ shelf_dist, northing
                      "log_wgt_total ~ s(shelf_distance_m, northing_m, k = 25) +
                                       s(t250m_Ano_fill, k = 5) +
                                       s(year_fac, bs = 're')",
                       
                      # # latxlon, S.V. t100, shelf dist
                      # "log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                      #                  s(easting_m, northing_m, by = t100m_Ano_fill, k = 10) +
                      #                  s(shelf_distance_m, k = 4) +
                      #                  s(year_fac, bs = 're')",

                      # S.V. T250 w/ easting/northing, shelf dist
                      "log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                                       s(easting_m, northing_m, by = t250m_Ano_fill, k = 10) +
                                       s(shelf_distance_m, k = 4) +
                                       s(year_fac, bs = 're')",
                      
                      # S.V. T250 w/ shelf_dist/northing
                      "log_wgt_total ~ s(shelf_distance_m, northing_m, k = 25) +
                                       s(shelf_distance_m, northing_m, by = t250m_Ano_fill, k = 10) +
                                       s(year_fac, bs = 're')",
                      
                      # Shelf Varying t250, shelf_dist, easting/northing - review question #2
                      "log_wgt_total ~ s(easting_m, northing_m, k = 25) +
                                       s(shelf_distance_m, k = 4) +
                                       s(shelf_distance_m, by = t250m_Ano_fill, k = 10) +
                                       s(year_fac, bs = 're')",
                      
                      # Shelf Varying t250, shelf_dist/northing
                      "log_wgt_total ~ s(shelf_distance_m, northing_m, k = 25) +
                                       s(shelf_distance_m, by = t250m_Ano_fill, k = 10) +
                                       s(year_fac, bs = 're')",
                      
                      # S.V. T250 w/ shelf_dist, northing
                      "log_wgt_total ~ s(northing_m, k = 20) +
                                       s(shelf_distance_m, k = 4) +
                                       s(shelf_distance_m, northing_m, by = t250m_Ano_fill, k = 10) +
                                       s(year_fac, bs = 're')"
)



############################################################################
# Total
############################################################################

# make formulas using strings above
formulas.bmass.total <- list()

for(i in 1:length(formula_strings)){
  formulas.bmass.total[[i]] <- formula(formula_strings[i])
}

formulas.pa.total <- lapply(formulas.bmass.total, FUN = function(x) update(x, present_total ~ .))

############################################################################
# Age 2 
############################################################################

formulas.bmass.age2 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age2 ~ .))
formulas.pa.age2 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age2 ~ .))

############################################################################
# Age 3-4
############################################################################

formulas.bmass.age3_4 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age3_4 ~ .))
formulas.pa.age3_4 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age3_4 ~ .))


############################################################################
# Age 5-20 
############################################################################

formulas.bmass.age5_20 <- lapply(formulas.bmass.total, FUN = function(x) update(x, log_age5_20 ~ .))
formulas.pa.age5_20 <- lapply(formulas.pa.total, FUN = function(x) update(x, present_age5_20 ~ .))



# Save formulas
############################################################################

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


if(!dir.exists("./gam/"))
  {dir.create("./gam/")}
save(formulas, file = "./gam/formulas.RData")



# end of script


