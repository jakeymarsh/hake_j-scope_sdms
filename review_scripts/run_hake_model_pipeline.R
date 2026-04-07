########################################################
## Master Analysis Pipeline
## Author: Jake Marshall  
## Date: September 2024
## Purpose: Runs complete hake modeling pipeline
########################################################


# Set to folder that has data and scripts folder
setwd("C:/Users/jakem/OneDrive/Desktop/GitHub/hake_j-scope_sdms/")


# Run Scripts in sequence
cat("Starting analysis pipeline...\n")

cat("1. Loading functions...\n")
source("./scripts/functions.R", echo=FALSE) 
source("./scripts/my_functions.R", echo=FALSE) 

cat("2. Making model formulas...\n") 
source("./scripts/make_model_formulas.R", echo=FALSE)

cat("3. Selecting best models for each age class...\n")
source("./scripts/model_selection.R", echo=FALSE)

cat("4. Generating output data...\n")
source("./scripts/generate_output.R", echo=FALSE)

cat("5. Generating visualizations for publication...\n")
source("./scripts/publish.R", echo=FALSE)


cat("Analysis pipeline complete!\n")





