#' #############################################################################
#' Project: Pacific Hake Distribution and Biomass Forecasting
#' Script: Pipeline
#' Author: Jake Marshall
#' Date: January 2026
#' Description: Orchestrates the end-to-end workflow for hake distribution 
#'              forecasting. This script sources all sub-modules in sequence, 
#'              from environmental data loading to the generation of 
#'              publication-quality transboundary biomass maps.
#' #############################################################################


# ------------------------------------------------------------------------------
# 1. Workspace Configuration
# ------------------------------------------------------------------------------
# Set the root directory for the GitHub repository
setwd("C:/Users/jakem/OneDrive/Desktop/GitHub/hake_j-scope_sdms/")

# Start Execution Timer
start_time <- Sys.time()
cat("Starting analysis pipeline at:", as.character(start_time), "\n")

# ------------------------------------------------------------------------------
# 2. Execution Modules
# ------------------------------------------------------------------------------

# Module 1: Environment Setup
# Loads libraries (mgcv, sf, terra) and custom Delta Model functions.
cat("Step 1: Loading functions and environmental assets...\n")
source("./scripts/functions.R", echo = FALSE) 
source("./scripts/load.R", echo = FALSE)

# Module 2: Model Definition
# Generates 40 candidate GAM formulas across four life-history stages.
cat("Step 2: Defining candidate GAM formulas...\n") 
source("./scripts/make_model_formulas.R", echo = FALSE)

# Module 3: Model Selection & LFO-CV
# Executes Leave-Future-Out cross-validation to select best biomass/presence models.
cat("Step 3: Performing LFO cross-validation and model selection...\n")
source("./scripts/model_selection.R", echo = FALSE)

# Module 4: Delta Model Simulations
# Runs 5,000 Monte Carlo simulations to generate 2024 distribution forecasts.
cat("Step 4: Generating spatially explicit Delta Model output...\n")
source("./scripts/generate_output.R", echo = FALSE)

# Module 5: Publication Assets
# Produces final heatmaps, partial effects plots, and EEZ aggregated signals.
cat("Step 5: Generating publication-quality visualizations...\n")
source("./scripts/publish.R", echo = FALSE)

# ------------------------------------------------------------------------------
# 3. Pipeline Completion
# ------------------------------------------------------------------------------
end_time <- Sys.time()
cat("Analysis pipeline complete!\n")
cat("Total Processing Time:", round(difftime(end_time, start_time, units="mins"), 2), "minutes\n")


# end of script

